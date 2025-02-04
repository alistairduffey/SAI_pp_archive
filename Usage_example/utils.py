## Utility funcitons for working with pp_archive data
import xarray as xr



def get_data(group, model, windows, table, variable, season='annual', mean_or_std='mean'):
    # group: 'ARISE', 'GeoMIP'
    
    # model: if group=='ARISE': UKESM1-0-LL or CESM2-WACCM; 
    ######## if group='GeoMIP': any of the 6 G6 models
    
    # windows: a list of window periods (and scenarios) over which to return values
    ######### if group =='ARISE': ['SSP245_background', 'SSP245_baseline', 'ARISE_assmt']
    ######### if group =='GeoMIP': ['G6sulfur_assmt', 'SSP245_baseline', 'SSP245_target', 'SSP585_background']
    
    # table: 'Amon', 'Omon', 'Lmon'
    # variable: many
    # season: 'annual', 'DJF', 'MAM', 'JJA', or 'SON'
    # mean_or_std: 'mean' or 'std'
    
    datas=[]
    for window in windows:
        path = '../pp_archive/{a}/{b}/maps/{c}/{d}/{e}/{f}/*_{g}_*.nc'.format(a=group, b=model, c=window,
                                                                           d=table, e=variable,
                                                                           f=mean_or_std, g=season)
        data = xr.open_mfdataset(path)
        datas.append(data)
    return datas



def CESMize_var_names(cmor_var):
    look_up = {'tas':'TREFHT',
               'pr':'PRECT'}
    
    return look_up[cmor_var]


def adjust_longitude(ds, lon_name='x'):
    """
    Adjust the longitude coordinate of an xarray dataset to -180 to 180.

    Parameters:
    - ds (xarray.Dataset): Input dataset with longitude coordinate.
    - lon_name (str): Name of the longitude coordinate in the dataset (default: 'x').

    Returns:
    - xarray.Dataset: Dataset with longitude adjusted to -180 to 180.
    """
    # Check if the longitude coordinate exists
    if lon_name not in ds.coords:
        raise ValueError(f"Longitude coordinate '{lon_name}' not found in the dataset.")
    
    # Extract longitude values
    lon = ds[lon_name]
    
    # Check if adjustment is necessary
    if lon.min() >= -180 and lon.max() <= 180:
        print("Longitude is already in the range -180 to 180. No adjustment needed.")
        return ds
    
    # Adjust longitude from 0 to 360 -> -180 to 180
    lon_new = ((lon + 180) % 360) - 180
    ds = ds.assign_coords({lon_name: lon_new})
    
    # Sort dataset by longitude for consistency
    ds = ds.sortby(lon_name)
    
    return ds