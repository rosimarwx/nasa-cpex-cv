import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import as_strided

def f_reflTexture(
    dbz: np.ndarray, 
    pix_rad: int, 
    dbz_base: float
) -> np.ndarray:
    """
    Calculates the texture of a 2D reflectivity field.

    This function is a Python translation of the f_reflTexture.m MATLAB
    function. It computes local variability in a sliding window after
    detrending the data.

    Args:
        dbz: 2D NumPy array of reflectivity data.
        pix_rad: The radius of the sliding window (window size = 2 * pix_rad + 1).
        dbz_base: A base reflectivity value to subtract from the data.

    Returns:
        A 2D NumPy array of the same shape as `dbz` containing the texture values.
    """
    if dbz.ndim != 2:
        raise ValueError("Input `dbz` must be a 2D array.")

    # Store original NaN locations to apply to the final output
    original_nan_mask = np.isnan(dbz)
    
    # 1. Pad the array with NaNs at the start and end of the second dimension
    # This handles the windowing at the edges of the data.
    dbz_padded = np.pad(
        dbz, 
        pad_width=((0, 0), (pix_rad, pix_rad)), 
        mode='constant', 
        constant_values=np.nan
    )

    # 2. Fill in missing data using linear interpolation.
    # Pandas is excellent for this, especially for handling start/end NaNs.
    df = pd.DataFrame(dbz_padded)
    # limit_direction='both' mimics MATLAB's 'EndValues','nearest'
    df_filled = df.interpolate(method='linear', axis=1, limit_direction='both')
    dbz_filled = df_filled.to_numpy()

    # 3. Adjust reflectivity with base value
    dbz_filled -= dbz_base
    
    # --- Vectorized Sliding Window Processing ---

    # 4. Create a view of all sliding windows using stride tricks
    # This avoids a slow Python loop and is highly memory-efficient.
    window_size = 2 * pix_rad + 1
    num_rows, padded_cols = dbz_filled.shape
    
    # The shape of the output will be (num_rows, num_original_cols, window_size)
    output_shape = (num_rows, dbz.shape[1], window_size)
    
    # Get the byte strides of the array
    row_stride, col_stride = dbz_filled.strides
    
    # The strides of the new view
    output_strides = (row_stride, col_stride, col_stride)
    
    dbz_windows = as_strided(dbz_filled, shape=output_shape, strides=output_strides)

    # 5. Calculate and remove the slope for all windows at once
    x1 = np.arange(window_size)
    
    # Use np.nansum for safe summation over the window axis (axis=2)
    sum_x = np.nansum(x1)
    sum_y = np.nansum(dbz_windows, axis=2)
    sum_xy = np.nansum(dbz_windows * x1, axis=2)
    sum_x2 = np.nansum(x1**2)
    
    N = np.sum(~np.isnan(dbz_windows), axis=2) # Count of non-NaN points in each window
    N[N == 0] = 1 # Avoid division by zero, result will be NaN anyway

    # Calculate slope (b) and intercept (a) for the linear fit in each window
    # The formulas are the standard solution for linear least squares
    denominator = (N * sum_x2 - sum_x**2).astype(float)
    # Avoid division by zero for windows with no variance in X (shouldn't happen here)
    denominator[denominator == 0] = np.nan 

    b = (N * sum_xy - sum_x * sum_y) / denominator
    a = (sum_y - b * sum_x) / N
    
    # Create the fitted line for each window using broadcasting
    # new_axis is needed to align arrays for element-wise operations
    fit_y = a[:, :, np.newaxis] + b[:, :, np.newaxis] * x1

    # 6. Remove the trend and add back the mean of the original window
    mean_y = np.nanmean(dbz_windows, axis=2)[:, :, np.newaxis]
    dbz_corr = dbz_windows - fit_y + mean_y
    dbz_corr[dbz_corr < 1] = 1 # Apply floor value

    # 7. Calculate texture for all windows
    # Texture = sqrt(std(detrended_data^2))
    variance = np.nanvar(dbz_corr**2, axis=2)
    tdbz = np.sqrt(variance)

    # 8. Create final output and apply original NaN mask
    dbz_text = tdbz
    dbz_text[original_nan_mask] = np.nan
    
    return dbz_text

