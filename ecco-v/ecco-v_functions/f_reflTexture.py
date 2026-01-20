import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
import matplotlib.pyplot as plt

def f_reflTexture(DBZ, pixRad, dbzBase):
    """
    Calculates reflectivity texture, translated directly from MATLAB.

    This version preserves the for-loop structure to ensure a one-to-one
    correspondence with the original MATLAB logic.

    Args:
        DBZ (np.ndarray): Input reflectivity data (2D array).
        pixRad (int): The pixel radius for the moving window.
        dbzBase (float): The base reflectivity value to subtract.

    Returns:
        np.ndarray: The calculated reflectivity texture.
    """
    # Initialize the output array with NaNs, same size as input
    dbzText = np.full(DBZ.shape, np.nan)

    # Pad data array at the start and end of the second dimension (columns)
    # This corresponds to MATLAB's padarray(DBZ, [0 pixRad], nan)
    dbzPadded = np.pad(
        DBZ, 
        pad_width=((0, 0), (pixRad, pixRad)), 
        mode='constant', 
        constant_values=np.nan
    )

#    # Fill in areas with no data (NaNs) using linear interpolation
#    # This is a direct equivalent of MATLAB's 
#    # fillmissing(..., 'linear', 2, 'EndValues', 'nearest')
#    # We process it row by row.
#    for i in range(dbzPadded.shape[0]):
#        s = pd.Series(dbzPadded[i, :])
#        # Interpolate fills NaNs in the middle
#        s.interpolate(method='linear', limit_direction='both', inplace=True)
#        # bfill/ffill handle any remaining NaNs at the very ends
#        s.bfill(inplace=True)
#        s.ffill(inplace=True)
#        dbzPadded[i, :] = s.values


# ... (previous code: dbzText initialization, padding, etc.) ...

# --- Corrected section for filling missing values ---
# This block correctly performs linear interpolation and extrapolation,
# matching MATLAB's fillmissing(...,'linear', 'EndValues', 'nearest').

    # We process it row by row.
    for i in range(dbzPadded.shape[0]):
        row_with_nans = dbzPadded[i, :]

        # Find the indices of the valid (non-NaN) data points
        all_indices = np.arange(len(row_with_nans))
        is_valid = ~np.isnan(row_with_nans)

        valid_indices = all_indices[is_valid]
        valid_values = row_with_nans[is_valid]

        # Handle edge cases where interpolation is not possible
        if len(valid_values) < 2:
            # If there are fewer than 2 data points, we can't create a line.
            # Fall back to filling with the single point, or leave as NaN if empty.
            if len(valid_values) == 1:
                dbzPadded[i, :] = valid_values[0] # Fill row with the single value
            continue # If row is all NaNs, leave it.

        # Create the interpolation function with extrapolation enabled
        f = interp1d(
            valid_indices,
            valid_values,
            kind='linear',            # Specifies linear interpolation
            bounds_error=False,       # Allows extrapolation
            fill_value="extrapolate"  # This is the key parameter
        )

        # Apply the function to all indices in the row
        dbzPadded[i, :] = f(all_indices)

# ... (rest of your code: subtracting dbzBase, the main loop, etc.) ...

    # Adjust reflectivity with base value
    dbzPadded = dbzPadded - dbzBase

    # Define the moving window width
    window_width = 2 * pixRad + 1
    
    # --- Loop through data points ---
    # Note: The loop iterates through each column of the *original* DBZ array
    # to calculate a texture value for each point.
    for ii in range(DBZ.shape[1]):
        # Pull out the correct window from the padded data
        dbzBlock = dbzPadded[:, ii : ii + window_width]

        # --- Calculate and remove the slope of reflectivity ---
        
        # Create the X matrix for the linear fit
        # Corresponds to: x1=1:size(dbzBlock,2); X=repmat(...)
        num_rows, num_cols = dbzBlock.shape
        x1 = np.arange(1, num_cols + 1)
        X = np.tile(x1, (num_rows, 1))
        
        # Calculate intermediate sums, ignoring NaNs
        # Corresponds to: sum(..., 2, 'omitnan')
        sumX = np.nansum(X, axis=1)
        sumY = np.nansum(dbzBlock, axis=1)
        sumXY = np.nansum(dbzBlock * X, axis=1)
        sumX2 = np.nansum(X**2, axis=1)
        sumY2 = np.nansum(dbzBlock**2, axis=1) # Calculated but not used, as in MATLAB

        # N is the number of non-NaN points in each row.
        # In this implementation, N should be constant (window_width) because
        # we have already filled the NaNs.
        N = num_cols

        # Calculate slope (b) and intercept (a) for each row
        denominator = (N * sumX2 - sumX**2)
        
        # Avoid division by zero
        # We need to cast the denominator to float to hold np.nan
        denominator = denominator.astype(float)
        denominator[denominator == 0] = np.nan

        a = (sumY * sumX2 - sumX * sumXY) / denominator
        b = (N * sumXY - sumX * sumY) / denominator

        # Recreate the linear fit line
        # We need to reshape a and b to (n, 1) to broadcast correctly with X (n, m)
        newY = a[:, np.newaxis] + b[:, np.newaxis] * X

        # Remove the slope (the linear trend)
        mean_dbzBlock = np.nanmean(dbzBlock, axis=1)
        dbzCorr = dbzBlock - newY + mean_dbzBlock[:, np.newaxis]
        
        # Set values less than 1 to 1
        dbzCorr[dbzCorr < 1] = 1

        # --- Calculate texture ---
        # NOTE: The original MATLAB code is sqrt(std(...)). This is unusual,
        # as std is already a square root. This calculates (variance)^0.25.
        # Translated literally as requested.
        # MATLAB's std() default is equivalent to numpy's std() with ddof=1.
        std_of_squares = np.nanstd(dbzCorr**2, axis=1, ddof=1)
        tdbz = np.sqrt(std_of_squares)

        # Assign the calculated texture vector to the output
        dbzText[:, ii] = tdbz

    # Final step: Ensure that any location that was originally NaN in the
    # input DBZ is also NaN in the output texture map.
    dbzText[np.isnan(DBZ)] = np.nan
    
    return dbzText
