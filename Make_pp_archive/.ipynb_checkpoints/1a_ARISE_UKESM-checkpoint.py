import os
import glob
import pandas as pd
import numpy as np
import xarray as xr
from xmip.preprocessing import rename_cmip6
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm
from utils import calc_weighted_spatial_means
import warnings
warnings.filterwarnings("ignore", module='xarray')

### options

# Model
model = 'UKESM1-0-LL'

# time-periods over which to take means
assessment_periods = {'SAI':slice('2050', '2069'),
                      'Background_warming':slice('2050', '2069'),
                      'Baseline':slice('2014', '2033'),
                      'Pre-industrial':slice('1850', '1900')}

# ARISE ensemble_members. We also only use these same members for the baseline
ens_mems = ['r1i1p1f2', 'r2i1p1f2', 'r3i1p1f2', 'r4i1p1f2', 'r8i1p1f2']


# all the CMIP7 BCVs which are monthly and on a single level:
# see markdown table below for details on variable meanings
vars_dict = {
             'prw':'Amon',
             'evspsbl':'Amon',
             'clivi':'Amon',
             'clt':'Amon',
             'clwvi':'Amon',
             'hfss':'Amon',
             'rlds':'Amon',
             'rsldscs':'Amon',
             'rlus':'Amon',
             'rlut':'Amon',
             'rlutcs':'Amon',
             'rsds':'Amon',
             'rsdscs':'Amon',
             'rsdt':'Amon',
             'rsus':'Amon',
             'rsuscs':'Amon',
             'rsut':'Amon',
             'rsutcs':'Amon',
             'pr':'Amon',
             'tas':'Amon',
             'uas':'Amon',
             'vas':'Amon',
             'hfls':'Amon',
             'hurs':'Amon',
             'huss':'Amon',
             'prc':'Amon',
             'prsn':'Amon',
             'ps':'Amon',
             'psl':'Amon',
             'sfcWind':'Amon',
             'tasmax':'Amon',
             'tasmin':'Amon',
             'tauu':'Amon',
             'tauv':'Amon',
             'ts':'Amon',
             'evspsblsoi':'Lmon',
             'lai':'Lmon',
             'mrfso':'Lmon',
             'mrro':'Lmon',
             'mrros':'Lmon',
             'mrso':'Lmon',
             'mrsos':'Lmon',
             'hfds':'Omon',
             'mlotst':'Omon',
             'sos':'Omon',
             'tauuo':'Omon',
             'tauvo':'Omon',
             'tos':'Omon',
             'zos':'Omon'
            }


# seasons
seasons = ['DJF', 'MAM', 'JJA', 'SON']

### get historical-SSP2-4.5 and combine along the time dimension

def get_historical_ssp245_ds(variable, table='Amon'):
    """ gets historical and SSP2-4.5 and combine along the time dimension
        returns dataset with 5 members, each running over the historical 
        1850-2014 and ssp245 2015-2100 """
    
    ds_list = []
    for es in ens_mems:
        path = '/badc/cmip6/data/CMIP6/ScenarioMIP/MOHC/UKESM1-0-LL/ssp245/{e}/{t}/{v}/*/latest/'.format(e=es, t=table,v=variable)
        ds = rename_cmip6(xr.open_mfdataset(path+'*.nc'))
        
        path_hist = glob.glob('/badc/cmip6/data/CMIP6/*/*/UKESM1-0-LL/historical/{e}/{t}/{v}/*/latest/'.format(
        t=table, v=variable, e=es))[0]
        ds_hist = rename_cmip6(xr.open_mfdataset(path_hist+'*.nc'))    
        ds = xr.concat([ds_hist, ds], dim='time')
        if 'height' in ds.variables:
            ds = ds.drop_vars('height')
        if 'type' in ds.variables:
            ds = ds.drop_vars('type')
        ds_list.append(ds)
    
    DS = xr.concat(ds_list, dim='Ensemble_member')
    return DS

## for ARISE
def get_ARISE_UKESM(variable='tas', table='Amon'):
    ds_list = []
    paths = glob.glob('/badc/deposited2022/arise/data/ARISE/MOHC/UKESM1-0-LL/arise-sai-1p5/*/{t}/{v}/*/*/'.format(
    t=table, v=variable))
    for path in paths:
        ds = rename_cmip6(xr.open_mfdataset(path+'*.nc'))
        if 'height' in ds.variables:
            ds = ds.drop_vars('height')
        if 'type' in ds.variables:
            ds = ds.drop_vars('type')
        ds_list.append(ds)
    DS = xr.concat(ds_list, dim='Ensemble_member')
    return DS

def get_time_period(ds, slice_label):
    ds_out = ds.sel(time=assessment_periods[slice_label])
    ds_out.attrs['t_bnds'] = str(assessment_periods[slice_label].start+'_'+assessment_periods[slice_label].stop)
    return ds_out

def process_and_save_maps(ds, ds_seasonal, var, table, label, seasons=seasons):
    """ 
    Inputs
    ds: a time resolved, quarterly resampled, spatial dataset, with an ensemble_member dimension
    label: 'SSP245_baseline', 'ARISE', or 'SSP245_background'. Defines naming of outputs. 
    
    Function saves the mean and standard deviation across the whole time+ens_mems combined dimension
    """
    path = '../pp_archive/ARISE/UKESM1-0-LL/maps/{l}/{t}/{v}/'.format(l=label, t=table, v=var)
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

for var in tqdm(vars_dict.keys()):

    print(var)
    # get data
    ds_ssp = get_historical_ssp245_ds(var, table=vars_dict[var])
    ds_arise = get_ARISE_UKESM(var, table=vars_dict[var])

    
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
    
    
    ds_ssp_seasonal = ds_ssp.resample(time="QS-DEC").mean()
    ds_arise_seasonal = ds_arise.resample(time="QS-DEC").mean()

    # process into time slice means:
    ssp_baseline, ssp_baseline_seasonal = get_time_period(ds_ssp, 'Baseline'), get_time_period(ds_ssp_seasonal, 'Baseline')
    ssp_background, ssp_background_seasonal = get_time_period(ds_ssp, 'Background_warming'), get_time_period(ds_ssp_seasonal, 'Background_warming')
    arise_assmt, arise_assmt_seasonal = get_time_period(ds_arise, 'SAI'), get_time_period(ds_arise_seasonal, 'SAI')

    process_and_save_maps(ssp_baseline, ssp_baseline_seasonal, 
                 var=var, table=vars_dict[var], 
                 label='SSP245_baseline', seasons=seasons)
    process_and_save_maps(ssp_background, ssp_background_seasonal, 
                     var=var, table=vars_dict[var], 
                     label='SSP245_background', seasons=seasons)
    process_and_save_maps(arise_assmt, arise_assmt_seasonal, 
                     var=var, table=vars_dict[var], 
                     label='ARISE_assmt', seasons=seasons)
