import numpy as np
from scipy import ndimage
from typing import Optional

def f_classSub(
    class_in: np.ndarray,
    asl: np.ndarray,
    topo: np.ndarray,
    melt_orig: np.ndarray,
    temp: np.ndarray,
    elev: Optional[np.ndarray],
    first_row: int,
    surf_alt_lim: float
) -> np.ndarray:
    """
    Refines a primary weather classification into detailed sub-categories.

    This is a Python translation of the f_classSub.m MATLAB function. It
    analyzes the properties of convective and stratiform regions to provide
    a more granular classification.

    Args:
        class_in: 2D array with initial classification codes (1: strat, 2: mixed, 3: conv).
        asl: 2D array of altitudes above sea level (meters).
        topo: 2D array of topography height (meters).
        melt_orig: 2D array of original melting layer data.
        temp: 2D array of temperature data (Celsius).
        elev: 1D array of elevation data related to aircraft position, or None if not available.
        first_row: The row index in the arrays corresponding to the aircraft's nadir.
        surf_alt_lim: Scalar altitude limit for defining "near surface" (meters).

    Returns:
        A 2D array with the detailed sub-classification codes.

    Classification Codes:
        - 14: Stratiform Low
        - 16: Stratiform Mid
        - 18: Stratiform High
        - 25: Mixed
        - 30: Convective (special case, e.g., near aircraft)
        - 32: Convective Elevated
        - 34: Convective Shallow
        - 36: Convective Mid
        - 38: Convective Deep
    """
    # Initialize the output array with NaNs
    class_sub = np.full(class_in.shape, np.nan)

    # --- Data Preparation ---
    
    # Calculate height above ground level (AGL)
    dist_agl = asl - topo

    # Create modifiable copies of melt and temp arrays
    melt = np.copy(melt_orig)
    temp_mod = np.copy(temp)
    
    # Modify melt and temp based on AGL, similar to the original logic
    # This seems to be a data conditioning step
    melt[dist_agl < 2000] = 9
    melt[np.isnan(melt_orig)] = np.nan
    temp_mod[(dist_agl < 4000) & (temp < -25)] = -25

    # --- Convective Processing ---

    # Create a mask for convective areas (classIn == 3)
    conv_mask = (class_in == 3)
    
    # Find and label each distinct convective area.
    # The structure=np.ones((3,3)) makes it an 8-connectivity problem, same as bwconncomp default.
    labeled_array, num_features = ndimage.label(conv_mask, structure=np.ones((3,3)))
    
    # Loop through each identified convective area (feature)
    for i in range(1, num_features + 1):
        # Create a mask for the single convective area we are currently processing
        area_mask = (labeled_array == i)
        
        # --- "Near Plane" Check ---
        plane_pix = 0
        if elev is not None:
            # Get row and column indices for the current area
            rows, cols = np.where(area_mask)
            
            # Find pixels in the area that are in the aircraft's row
            is_plane_row_mask = (rows == first_row)
            plane_pix_total = np.sum(is_plane_row_mask)
            
            if plane_pix_total > 0:
                # Get AGL for pixels in the plane's row
                alt_diff_plane_row = dist_agl[area_mask][is_plane_row_mask]
                
                # Get column indices for these specific pixels
                cols_plane_row = cols[is_plane_row_mask]
                
                # Get elevation data for these pixels
                elev_plane_pix = elev[cols_plane_row]

                # Find pixels considered too low relative to the plane
                alt_low = np.sum((alt_diff_plane_row < 500) & (elev_plane_pix > 0))
                plane_pix = plane_pix_total - alt_low
        
        # --- "Near Surface" Check ---
        asl_area = dist_agl[area_mask]
        near_surf_pix = np.sum(asl_area < (500 + surf_alt_lim))
        
        # --- Apply Classification Logic for the Convective Cell ---
        if near_surf_pix == 0:  # Not near surface: Elevated Convection
            # if plane_pix > 10 and np.median(elev[cols]) > 0: # This check is complex to translate directly without full context of 'elev'
            # For now, implementing the core logic. The 'plane' override is specific.
            # Assuming 'cols' are from the entire area_mask for the median check.
            _, area_cols = np.where(area_mask)
            if elev is not None and plane_pix > 10 and np.median(elev[area_cols]) > 0:
                 class_sub[area_mask] = 30
            else:
                 class_sub[area_mask] = 32
        else:  # Near surface: Shallow, Mid, or Deep Convection
            melt_max = np.nanmax(melt[area_mask])
            if melt_max < 15:  # Below melting layer: Shallow
                if elev is not None and plane_pix > 10:
                    class_sub[area_mask] = 30
                else:
                    class_sub[area_mask] = 34
            else: # Above melting layer
                min_temp = np.nanmin(temp_mod[area_mask])
                if min_temp >= -25:  # Below divergence level: Mid
                    if elev is not None and plane_pix > 10:
                        class_sub[area_mask] = 30
                    else:
                        class_sub[area_mask] = 36
                else:  # Deep
                    _, area_cols = np.where(area_mask)
                    if elev is not None and plane_pix > 10 and np.median(elev[area_cols]) > 0:
                        class_sub[area_mask] = 30
                    else:
                        class_sub[area_mask] = 38

    # --- Stratiform and Mixed Processing (using boolean masking) ---
    
    # Stratiform Low
    strat_low_mask = (class_in == 1) & (melt < 15)
    class_sub[strat_low_mask] = 14
    
    # Stratiform Mid
    strat_mid_mask = (class_in == 1) & (melt >= 15) & (temp_mod >= -25)
    class_sub[strat_mid_mask] = 16
    
    # Stratiform High
    strat_high_mask = (class_in == 1) & (melt >= 15) & (temp_mod < -25)
    class_sub[strat_high_mask] = 18

    # Mixed
    mixed_mask = (class_in == 2)
    class_sub[mixed_mask] = 25
    
    return class_sub
