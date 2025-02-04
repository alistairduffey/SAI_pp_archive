### functions to account for month length variation in means:
import xarray as xr
import numpy as np

def season_mean(ds):
    ### this function is from Joe Hamman's xarray docs 
    ### https://docs.xarray.dev/en/stable/examples/monthly-means.html
    
    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = ds.time.dt.days_in_month

    # Calculate the weights by grouping by 'time.season'
    weights = (
        month_length.groupby("time.season") / month_length.groupby("time.season").sum()
    )

    # Test that the sum of the weights for each season is 1.0
    np.testing.assert_allclose(weights.groupby("time.season").sum().values, np.ones(4))

    # Calculate the weighted average
    return (ds * weights).groupby("time.season").sum(dim="time")


def weighted_seasonal_resample(ds):
    """
    weight by days in each month
    adapted from NCAR docs 
    https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/
    """
    # Make a DataArray with the number of days in each month, size = len(time)
    month_length = ds.time.dt.days_in_month

    # Calculate the weights by grouping by 'time.season'
    wgts = (
        month_length.groupby("time.season") / month_length.groupby("time.season").sum()
    )

    # Test that the sum of the weights for each season is 1.0
    np.testing.assert_allclose(wgts.groupby("time.season").sum().values, np.ones(4))
    
    numerator = (ds * wgts).resample(time="QS-DEC").sum(dim="time")
    denominator = wgts.resample(time="QS-DEC").sum(dim="time")

    return numerator/denominator

def weighted_annual_resample(ds):
    """
    weight by days in each month
    adapted from NCAR docs 
    https://ncar.github.io/esds/posts/2021/yearly-averages-xarray/
    """
    # Determine the month length
    month_length = ds.time.dt.days_in_month

    # Calculate the weights
    wgts = month_length.groupby("time.year") / month_length.groupby("time.year").sum()

    # Make sure the weights in each year add up to 1
    np.testing.assert_allclose(wgts.groupby("time.year").sum(xr.ALL_DIMS), 1.0)

    numerator = (ds * wgts).resample(time="YS").sum(dim="time")
    denominator = wgts.resample(time="YS").sum(dim="time")

    return numerator/denominator


def calc_weighted_spatial_means(ds, var_name, land_mask=None):
    """
    Get the area-weighted and land-weighted mean of a spatial xarray dataset.

    Parameters:
    -----------
    ds : xarray.Dataset or xarray.DataArray
        The dataset containing the variable of interest.
    var_name : str
        Name of the variable for which to compute weighted means.
    land_mask : xarray.DataArray, optional
        A land frac dataset (100 for fully land, 0 for ocean) with dimensions (y, x).
    
    Returns:
    --------
    dict
        Dictionary with area-weighted and land-weighted means.
    """

    # Extract the variable of interest
    da = ds[var_name] if isinstance(ds, xr.Dataset) else ds

    # use cosine of latitude for area weighting
    lat_radians = np.deg2rad(da.y)
    lat_weights = np.cos(lat_radians)
    
    # Compute area-weighted mean efficiently using xarray's `weighted()`
    area_mean = da.weighted(lat_weights).mean(dim=("y", "x")).to_dataset(name=var_name)

    # Compute land-weighted mean (if land mask is provided)
    land_mean, ocean_mean = None, None
    if land_mask is not None:
        land_weighted = lat_weights * land_mask  # Keeps lat_weights as 1D
        ocean_weighted = lat_weights * (100 - land_mask)

        land_mean = da.weighted(land_weighted).mean(dim=("y", "x")).to_dataset(name=var_name)
        ocean_mean = da.weighted(ocean_weighted).mean(dim=("y", "x")).to_dataset(name=var_name)

    return {
            "area_weighted_mean": area_mean,
            "land_weighted_mean": land_mean if land_mask is not None else None,
            "ocean_weighted_mean": ocean_mean if land_mask is not None else None,
            }