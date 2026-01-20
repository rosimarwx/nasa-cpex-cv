import xarray as xr
import dask.array as da
import numpy as np

def calculate_texture_chunk(dbz_padded_chunk, pixRad, dbzBase, window_width):
    """
    Calculates reflectivity texture on a single chunk of data using a vectorized
    (non-looping) approach.
    """
    # Guard clause for small chunks remains important
    output_width = dbz_padded_chunk.shape[1] - (window_width - 1)
    if output_width <= 0:
        return np.empty((dbz_padded_chunk.shape[0], 0), dtype=np.float64)

    # Adjust reflectivity with base value
    dbz_padded_chunk = dbz_padded_chunk - dbzBase
    
    # --- Vectorization using stride tricks ---
    # Create a view of all sliding windows. This is the key step that
    # replaces the for-loop. The shape will be (num_rows, num_windows, window_width)
    shape = (dbz_padded_chunk.shape[0], output_width, window_width)
    strides = (dbz_padded_chunk.strides[0], dbz_padded_chunk.strides[1], dbz_padded_chunk.strides[1])
    dbzBlock = np.lib.stride_tricks.as_strided(dbz_padded_chunk, shape=shape, strides=strides)

    # Now, all operations are performed on the window axis (axis=2)
    num_rows, num_cols = dbzBlock.shape[0], dbzBlock.shape[2] # num_rows=1, num_cols=11
    x1 = np.arange(1, num_cols + 1)
    # Reshape X for broadcasting with the (num_rows, num_windows, window_width) dbzBlock
    X = x1.reshape(1, 1, -1)

    sumX = np.nansum(X, axis=2)
    sumY = np.nansum(dbzBlock, axis=2)
    sumXY = np.nansum(dbzBlock * X, axis=2)
    sumX2 = np.nansum(X**2, axis=2)
    
    N = np.sum(~np.isnan(dbzBlock), axis=2)
    denominator = (N * sumX2 - sumX**2)
    denominator[denominator == 0] = np.nan

    a = (sumY * sumX2 - sumX * sumXY) / denominator
    b = (N * sumXY - sumX * sumY) / denominator

    # Add a new axis for broadcasting against the window dimension
    newY = a[..., np.newaxis] + b[..., np.newaxis] * X
    mean_dbzBlock = np.nanmean(dbzBlock, axis=2)
    dbzCorr = dbzBlock - newY + mean_dbzBlock[..., np.newaxis]
    dbzCorr[dbzCorr < 1] = 1

    std_of_squares = np.nanstd(dbzCorr**2, axis=2, ddof=1)
    tdbz = np.sqrt(std_of_squares)
    
    return tdbz # The result is already the correct shape

# This is the main function you will call with your Dask array.
def f_reflTexture_dask(dbz_da, pixRad, dbzBase):
    """
    Calculates reflectivity texture on a Dask-backed xarray.DataArray.

    Args:
        dbz_da (xr.DataArray): Input reflectivity data, likely with dims ('vertical', 'horizontal').
                               The moving window is applied along the horizontal dimension.
        pixRad (int): The pixel radius for the moving window.
        dbzBase (float): The base reflectivity value to subtract.

    Returns:
        xr.DataArray: A Dask-backed DataArray of the calculated reflectivity texture.
    """
    # Assume the horizontal dimension is the last one.
    ##horizontal_dim = dbz_da.dims[-1]
    
    # RRB: I believe step 1 isn't needed for model data
    # --- Step 1: Fill missing values using Xarray's Dask-aware method ---
    # This replaces the first for-loop in your original function.
    # print("Step 1: Interpolating NaN values...")
    # dbz_filled = dbz_da.interpolate_na(
    #    dim=horizontal_dim, 
    #    method="linear",
    #    kwargs={"fill_value": "extrapolate"}
    #)
    dbz_filled = dbz_da

    # --- Step 2: Apply the texture calculation using Dask's map_overlap ---
    window_width = 2 * pixRad + 1
    
    # The output at column `i` depends on inputs from `i` to `i + window_width - 1`.
    # This means we need an overlap on the right side of each chunk.
    # The depth is the number of extra cells needed from the neighbor.
    overlap_depth = window_width - 1

    # Apply the custom chunk-wise function.
    # The result is a raw Dask array, without coordinates.
    dbz_texture_dask_array = da.map_overlap(
        calculate_texture_chunk,
        dbz_filled.data,  # Pass the underlying Dask array
        depth={1: overlap_depth},  # Depth on axis 1 (horizontal), 0 on axis 0
        boundary='none', # Use 'none' as we handle padding manually with NaNs if needed at edges
        trim=True, # Trim the overlap from the output of each chunk
        dtype=dbz_filled.dtype,
        pixRad=pixRad, # Additional keyword arguments for the function
        dbzBase=dbzBase,
        window_width=window_width
    )
    
    # --- Step 3: Wrap the result back into an xarray.DataArray and mask it ---

    computed_coords = {
        dim: coord.compute() 
        for dim, coord in dbz_da.coords.items() 
        if dim in dbz_da.dims
    }
    
    # Create the final DataArray, preserving original coordinates and dimensions.
    dbz_texture_da = dbz_texture_dask_array
#    dbz_texture_da = xr.DataArray(
#        dbz_texture_dask_array,
#        coords = computed_coords,
        # dims=dbz_da.dims,
#        name='reflectivity_texture'
#    )
    
    # Final step: Ensure that any location that was originally NaN in the
    # input is also NaN in the output.
    return dbz_texture_da#.where(~dbz_da.isnull())

