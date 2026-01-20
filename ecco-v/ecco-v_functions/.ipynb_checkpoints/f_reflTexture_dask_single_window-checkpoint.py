import xarray as xr
import dask.array as da
import numpy as np

def calculate_texture_single_window(window, dbzBase):
    """
    Calculates texture on a single 1D window of data.
    This is the core function to be used with apply_ufunc.
    """
    # Adjust reflectivity with base value
    dbzBlock = window - dbzBase
    
    # --- The same linear fit logic, but on a 1D array ---
    num_cols = len(dbzBlock)
    x1 = np.arange(1, num_cols + 1)
    
    sumX = np.nansum(x1)
    sumY = np.nansum(dbzBlock)
    sumXY = np.nansum(dbzBlock * x1)
    sumX2 = np.nansum(x1**2)
    
    N = np.sum(~np.isnan(dbzBlock))
    if N < 2: # Cannot do a linear fit on less than 2 points
        return np.nan

    denominator = (N * sumX2 - sumX**2)
    if denominator == 0:
        return np.nan

    a = (sumY * sumX2 - sumX * sumXY) / denominator
    b = (N * sumXY - sumX * sumY) / denominator

    newY = a + b * x1
    mean_dbzBlock = np.nanmean(dbzBlock)
    dbzCorr = dbzBlock - newY + mean_dbzBlock
    dbzCorr[dbzCorr < 1] = 1

    std_of_squares = np.nanstd(dbzCorr**2, ddof=1)
    tdbz = np.sqrt(std_of_squares)
    
    return tdbz

def f_reflTexture_dask(dbz_da, pixRad, dbzBase):
    """
    Calculates reflectivity texture using the robust xarray.rolling().construct() method.
    """
    window_width = 2 * pixRad + 1
    horizontal_dim = dbz_da.dims[-1]

    # Step 1: Pad the array with NaNs on the right side.
    # This ensures the window at each point `i` includes data from `i` to `i + window_width -1`,
    # matching the original logic.
    pad_width = window_width - 1
    padded_da = dbz_da.pad(pad_width={horizontal_dim: (0, pad_width)}, constant_values=np.nan)

    # Step 2: Create a DataArray where each element contains the entire rolling window.
    # The `center=False` makes the window start at each coordinate.
    # Create a dictionary where the key is the *value* of horizontal_dim ('nCells')
    window_dim_dict = {horizontal_dim: "window_dim"}
    
    # Unpack the dictionary using ** to pass `nCells="window_dim"` as the keyword argument
    rolling_windows = padded_da.rolling(
        {horizontal_dim: window_width}, center=False
    ).construct(**window_dim_dict)
    # rolling_windows = padded_da.rolling({horizontal_dim: window_width}, center=False).construct(horizontal_dim="window_dim")
    
    # Step 3: Apply our simple window function across all windows using apply_ufunc.
    # This is a highly optimized and stable Dask operation.
    reflectivity_texture = xr.apply_ufunc(
        calculate_texture_single_window,
        rolling_windows,
        input_core_dims=[["window_dim"]], # The function operates on the "window_dim"
        dask="parallelized",
        output_dtypes=[dbz_da.dtype],
        kwargs={'dbzBase': dbzBase} # Pass dbzBase as a keyword argument
    )

    # The output of rolling is smaller, so we need to align it back to the original coordinates
    # This ensures the output has the exact same shape and dims as the input dbz_da
    # aligned_texture = reflectivity_texture.reindex_like(dbz_da)
        # Instead of reindex_like, we slice the result back to the original input size,
    # effectively trimming off the padding.
    slicers = {dim: slice(None, size) for dim, size in dbz_da.sizes.items()}
    aligned_texture = reflectivity_texture.isel(**slicers)

    return aligned_texture.where(~dbz_da.isnull())