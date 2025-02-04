import os
import glob
import pandas as pd
import numpy as np
import xarray as xr
from xmip.preprocessing import rename_cmip6
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm
from utils import weighted_seasonal_resample, weighted_annual_resample, calc_weighted_spatial_means
import warnings
warnings.filterwarnings("ignore", module='xarray')

### options

# Model
model = 'CESM2-WACCM'

# time-periods over which to take means
assessment_periods = {'SAI':slice('2050', '2069'),
                      'Background_warming':slice('2050', '2069'),
                      'Baseline':slice('2020', '2039'), # see 
                      'Pre-industrial':slice('1850', '1900')}

## note that CESM ARISE data is not CMOR-ized
cesm_vars = ['TREFHT', 'PRECT', 'PSL', 
             'LHFLX', 'TS', 'FSDS', 
             'TREFHTMX', 'TREFHTMN', 'RHREFHT', 
             'AICE']


ensemble_members_CESM = ['001', '002', '003', '004', '005',
                         '006', '007', '008', '009', '010']

experiment_paths_CESM = {'ARISE':'/gws/nopw/j04/cpom/aduffey/ARISE/CESM/arise/',
                         'SSP245':'/gws/nopw/j04/cpom/aduffey/ARISE/CESM/ssp245/'}

# seasons
seasons = ['DJF', 'MAM', 'JJA', 'SON']


### get data

def get_ssp245_or_arise_ds(variable, scenario, ens_mems):
    """ returns dataset with 10 members, each running over ssp245 2015-2100 """
    
    ds_list = []
    for es in ens_mems:
        path = experiment_paths_CESM[scenario] + '{v}/'.format(v=variable)
        ds = rename_cmip6(xr.open_mfdataset(path+'*.{}.c*.nc'.format(es)))
        if 'height' in ds.variables:
            ds = ds.drop_vars('height')
        if 'type' in ds.variables:
            ds = ds.drop_vars('type')
            
        ## rename 
        ds_list.append(ds)
    
    DS = xr.concat(ds_list, dim='Ensemble_member')
    return DS


def get_time_period(ds, slice_label):
    ds_out = ds.sel(time=assessment_periods[slice_label])
    ds_out.attrs['t_bnds'] = str(assessment_periods[slice_label].start+'_'+assessment_periods[slice_label].stop)
    return ds_out

def process_and_save_maps(ds, ds_seasonal, var, table, label, seasons):
    """ 
    Inputs
    ds: a time resolved, quarterly resampled, spatial dataset, with an ensemble_member dimension
    label: 'SSP245_baseline', 'ARISE', or 'SSP245_background'. Defines naming of outputs. 
    
    Function saves the mean and standard deviation across the whole time+ens_mems combined dimension
    """
    path = '../pp_archive/ARISE/CESM2-WACCM/maps/{l}/{t}/{v}/'.format(l=label, t=table, v=var)
    os.makedirs(path+'/std/', exist_ok=True)
    os.makedirs(path+'/mean/', exist_ok=True)
    
    t_bnds = ds_seasonal.t_bnds
    for season in seasons:
        ds_season = ds_seasonal.where(ds_seasonal.time.dt.season == season, drop=True)
        std = ds_season.std(dim=['time', 'Ensemble_member'])
        mean = ds_season.mean(dim=['time', 'Ensemble_member'])
        
        std.to_netcdf(path + '/std/' +'{v}_{l}_{s}_{t}_std.nc'.format(v=var, l=label, s=season, t=t_bnds))
        mean.to_netcdf(path + '/mean/' + '{v}_{l}_{s}_{t}_mean.nc'.format(v=var, l=label, s=season, t=t_bnds))

    t_bnds = ds.t_bnds
    # repeat for the annual mean:
    std = ds.std(dim=['time', 'Ensemble_member'])
    mean = ds.mean(dim=['time', 'Ensemble_member'])
    
    std.to_netcdf(path + '/std/' + '{v}_{l}_annual_{t}_std.nc'.format(v=var, l=label, t=t_bnds))
    mean.to_netcdf(path + '/mean/' + '{v}_{l}_annual_{t}_mean.nc'.format(v=var, l=label, t=t_bnds))
    return

def save_timeseries(ds, var, table, label, scenario):
    """ 
    Inputs
    ds: a monthly, non-spatial dataset
    label: 'land', 'ocean', 'global'
    
    Function saves the timeseries for each ens_mems
    """
    
    outpath = '../pp_archive/ARISE/{m}/timeseries/{l}/{t}/{v}/'.format(m=model, l=label, t=table, v=var)
    os.makedirs(outpath, exist_ok=True)
    ds.to_netcdf(outpath + '{v}_{l}_{s}_timeseries.nc'.format(v=var, l=label, s=scenario))
    return


# MAIN

for var in tqdm(cesm_vars):

    print(var)
    # get data
    ds_ssp = get_ssp245_or_arise_ds(var, 'SSP245', ensemble_members_CESM)
    ds_arise = get_ssp245_or_arise_ds(var, 'ARISE', ensemble_members_CESM)
    
    ## first save the global, land and ocean mean time series
    land_fraction = xr.open_dataset('../pp_archive/fx/{m}/sftlf/{m}_sftlf.nc'.format(m=model)).sftlf
    i=0
    scens = ['SSP245', 'ARISE']
    for run_ds in [ds_ssp.copy(), ds_arise.copy()]:
        scen = scens[i]
        means = calc_weighted_spatial_means(run_ds, var_name=var, land_mask=land_fraction)
        save_timeseries(means["land_weighted_mean"], var=var, scenario=scen, table='Amon', label='land')
        save_timeseries(means["ocean_weighted_mean"], var=var, scenario=scen, table='Amon', label='ocean')
        save_timeseries(means["area_weighted_mean"], var=var, scenario=scen, table='Amon', label='global')
        i=i+1

    

    ### resample to annual resolution, weighting by month length
    ds_ssp_annual = weighted_annual_resample(ds_ssp)
    ds_arise_annual = weighted_annual_resample(ds_arise)

    ### resample to seasonal resolution, weighting by month length
    ds_ssp_seasonal = weighted_seasonal_resample(ds_ssp)
    ds_arise_seasonal = weighted_seasonal_resample(ds_arise)

    # process into time slice means:
    ssp_baseline, ssp_baseline_seasonal = get_time_period(ds_ssp, 'Baseline'), get_time_period(ds_ssp_seasonal, 'Baseline')
    ssp_background, ssp_background_seasonal = get_time_period(ds_ssp, 'Background_warming'), get_time_period(ds_ssp_seasonal, 'Background_warming')
    arise_assmt, arise_assmt_seasonal = get_time_period(ds_arise, 'SAI'), get_time_period(ds_arise_seasonal, 'SAI')

    process_and_save_maps(ssp_baseline, ssp_baseline_seasonal, 
                          var=var, table='Amon', 
                          label='SSP245_baseline', seasons=seasons)
    process_and_save_maps(ssp_background, ssp_background_seasonal, 
                          var=var, table='Amon', 
                          label='SSP245_background', seasons=seasons)
    process_and_save_maps(arise_assmt, arise_assmt_seasonal, 
                          var=var, table='Amon', 
                          label='ARISE_assmt', seasons=seasons)
