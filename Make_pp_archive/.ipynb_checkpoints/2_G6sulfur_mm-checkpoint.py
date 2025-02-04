### postprocess GeoMIP

#Alistair Duffey, october 2024
#All models with G6sulfur runs
#N.B.: runs over the CEDA archive data structure 
#N.B: MIP-ESM1-2-HR: no SSP245/585 runs of the r3 ensemble member, so only keep two. 


#### output structure 
## pp_archive/GeoMIP/model/scenario/ensemble_member/maps/table/variable/file.nc
## note that while the number of members varies between 1 and a few, we save as ens mean and std in each case, but also save the number of members. 

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

models = ['IPSL-CM6A-LR', 'UKESM1-0-LL', 'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'CESM2-WACCM', 'CNRM-ESM2-1']

scenarios = ['G6sulfur', 'ssp245', 'ssp585']

scenario_types = {'G6sulfur':'GeoMIP',
                  'ssp245':'ScenarioMIP',
                  'ssp585':'ScenarioMIP'}


# time-periods over which to take means
assessment_periods = {'Future':slice('2080', '2099'),
                      'Baseline':slice('2015', '2034')}

## we don't use ensemble member inputs here - instead search over available members in each case

# note that we only use a member when that same label is also there for the ssps - this isn't ideal. #To improve

ensemble_members = {
                    'IPSL-CM6A-LR':['r1i1p1f1'],
                    'UKESM1-0-LL':['r1i1p1f2', 'r4i1p1f2', 'r8i1p1f2'],
                    'MPI-ESM1-2-HR':['r1i1p1f1', 'r2i1p1f1'],#, 'r3i1p1f1'],
                    'MPI-ESM1-2-LR':['r1i1p1f1', 'r2i1p1f1', 'r3i1p1f1'],
                    'CESM2-WACCM':['r1i1p1f2', 'r2i1p1f2'],
                    'CNRM-ESM2-1':['r1i1p1f2', 'r2i1p1f2', 'r3i1p1f2']
                    }


## for cesm2, exclude some G6 files, because the archive has overlapping time data
excluded_files = ['tas_Amon_CESM2-WACCM_G6sulfur_r1i1p1f2_gn_202001-206912.nc',
                  'tas_Amon_CESM2-WACCM_G6sulfur_r2i1p1f2_gn_201501-206412.nc',
                  'tas_Amon_CESM2-WACCM_G6sulfur_r2i1p1f2_gn_206501-210012.nc',
                  'tas_Amon_CESM2-WACCM_G6sulfur_r1i1p1f2_gn_207001-210012.nc']

# top 10 CMIP6 most downloaded variables:
# see markdown table below for details on variable meanings
vars_dict = {
             'evspsbl':'Amon',
             'rsds':'Amon',
             'pr':'Amon',
             'tas':'Amon',
             'hurs':'Amon',
             'psl':'Amon',
             'tasmax':'Amon',
             'tasmin':'Amon',
             'ts':'Amon',
            'siconca':'SImon',
            }


# seasons
seasons = ['DJF', 'MAM', 'JJA', 'SON']

### get data
def get_data_ssp_G6(model, scenario, variable, table):
    root = '/badc/cmip6/data/CMIP6/' # CEDA archive root
    
    ens_mems_mod = ensemble_members[model]
    ds_list = []
    
    for es in ens_mems_mod:
        
        ##############################################################
        ## following lines are a hack to get around the fact that CESM2 
        ## has f1 variants for ssps, and f2 variants for geoMIP
        if model == 'CESM2-WACCM':
            if scenario != 'G6sulfur':
                es = es.replace('p1f2', 'p1f1')
        ############################################################
        directory = root + '*/*/{m}/{s}/{e}/{t}/{v}/*/latest/'.format(m=model, s=scenario, t=table, v=variable, e=es)
        try: 
            path = glob.glob(directory)[0]
            ds = rename_cmip6(xr.open_mfdataset(path+'*.nc', use_cftime=True))
        
        except:
            print('main route failed - trying CESM workaround')
            if model == 'CESM2-WACCM':
                if scenario == 'G6sulfur':
                    files = []
                    for x in os.listdir(glob.glob(directory)[0]):
                        files.append(glob.glob(directory)[0] + x)
                    # drop extraneous data from cesm G6sulfur runs:
                    for file in excluded_files:
                        files = [s for s in files if not file in s]
                    ds = rename_cmip6(xr.open_mfdataset(files, use_cftime=True))
            
        if 'height' in ds.variables:
            ds = ds.drop_vars('height')
        if 'type' in ds.variables:
            ds = ds.drop_vars('type')
            
        ## rename 
        ds_list.append(ds)
    
    DS = xr.concat(ds_list, dim='Ensemble_member')
    DS = DS.assign_coords({'Ensemble_member':ens_mems_mod})
    return DS

def get_time_period(ds, slice_label):
    ds_out = ds.sel(time=assessment_periods[slice_label])
    ds_out.attrs['t_bnds'] = str(assessment_periods[slice_label].start+'_'+assessment_periods[slice_label].stop)
    return ds_out

def process_and_save_maps(ds, ds_seasonal, 
                          var, table, 
                          label, seasons=seasons):
    """ 
    Inputs
    ds: a time resolved, quarterly resampled, spatial dataset, with an ensemble_member dimension
    label: 'SSP245_baseline', 'G6sulfur', 'SSP245_target', 'SSP585_target'. Defines naming of outputs. 
    
    Function saves the mean and time standard deviation for each ens_mems
    """
    
    outpath = '../pp_archive/GeoMIP/{m}/maps/{l}/{t}/{v}/'.format(m=model, l=label, 
                                                               t=table, v=var)
    
    os.makedirs(outpath+'/std/', exist_ok=True)
    os.makedirs(outpath+'/mean/', exist_ok=True)
    
    t_bnds = ds_seasonal.t_bnds
    for season in seasons:
        ds_season = ds_seasonal.where(ds_seasonal.time.dt.season == season, drop=True)
        ds_season['N_members'] = len(ds_season.Ensemble_member.values)
        std = ds_season.std(dim=['time', 'Ensemble_member'])
        mean = ds_season.mean(dim=['time', 'Ensemble_member'])
        
        std.to_netcdf(outpath + '/std/' +'{v}_{l}_{s}_{t}_std.nc'.format(v=var, l=label,
                                                                           s=season, t=t_bnds))
        mean.to_netcdf(outpath + '/mean/' + '{v}_{l}_{s}_{t}_mean.nc'.format(v=var, l=label, 
                                                                               s=season, t=t_bnds))

    t_bnds = ds.t_bnds
    # repeat for the annual mean:
    ds['N_members'] = len(ds.Ensemble_member.values)
    std = ds.std(dim=['time', 'Ensemble_member'])
    mean = ds.mean(dim=['time', 'Ensemble_member'])
    
    std.to_netcdf(outpath + '/std/' + '{v}_{l}_annual_{t}_std.nc'.format(v=var, l=label, 
                                                                           t=t_bnds))
    mean.to_netcdf(outpath + '/mean/' + '{v}_{l}_annual_{t}_mean.nc'.format(v=var, l=label, 
                                                                             t=t_bnds))
    return

def save_timeseries(ds, var, table, label, scenario):
    """ 
    Inputs
    ds: a monthly, non-spatial dataset
    label: 'land', 'ocean', 'global'
    
    Function saves the timeseries for each ens_mems
    """
    
    outpath = '../pp_archive/GeoMIP/{m}/timeseries/{l}/{t}/{v}/'.format(m=model, l=label, t=table, v=var)
    os.makedirs(outpath, exist_ok=True)
    ds.to_netcdf(outpath + '{v}_{l}_{s}_timeseries.nc'.format(v=var, l=label, s=scenario))
    return

# MAIN
for model in models:
    print(model)
    for var in tqdm(vars_dict.keys()):
    
        print(var)
        # get data
        try:
            ds_ssp245 = get_data_ssp_G6(model=model, scenario='ssp245',
                                        variable=var, table=vars_dict[var])
            ds_ssp585 = get_data_ssp_G6(model=model, scenario='ssp585',
                                        variable=var, table=vars_dict[var])
            ds_G6sulfur = get_data_ssp_G6(model=model, scenario='G6sulfur',
                                        variable=var, table=vars_dict[var])

            
            ## first save the global, land and ocean mean time series
            land_fraction = xr.open_dataset('../pp_archive/fx/{m}/sftlf/{m}_sftlf.nc'.format(m=model)).sftlf
            i=0
            scens = ['SSP245', 'SSP585', 'G6sulfur']
            for run_ds in [ds_ssp245.copy(), ds_ssp585.copy(), ds_G6sulfur.copy()]:
                scen = scens[i]
                means = calc_weighted_spatial_means(run_ds, var_name=var, land_mask=land_fraction)
                save_timeseries(means["land_weighted_mean"], var=var, scenario=scen, table='Amon', label='land')
                save_timeseries(means["ocean_weighted_mean"], var=var, scenario=scen, table='Amon', label='ocean')
                save_timeseries(means["area_weighted_mean"], var=var, scenario=scen, table='Amon', label='global')
                i=i+1

            
            ### resample to annual resolution, weighting by month length
            ds_ssp245_annual = weighted_annual_resample(ds_ssp245)
            ds_ssp585_annual = weighted_annual_resample(ds_ssp585)
            ds_G6sulfur_annual = weighted_annual_resample(ds_G6sulfur)
            
            ### resample to seasonal resolution, weighting by month length
            ds_ssp245_seasonal = weighted_seasonal_resample(ds_ssp245)
            ds_ssp585_seasonal = weighted_seasonal_resample(ds_ssp585)
            ds_G6sulfur_seasonal = weighted_seasonal_resample(ds_G6sulfur)
        
            # process into time slice means:
            ssp245_baseline, ssp245_baseline_seasonal = get_time_period(ds_ssp245_annual, 'Baseline'), get_time_period(ds_ssp245_seasonal, 'Baseline')
            ssp585_background, ssp585_background_seasonal = get_time_period(ds_ssp585_annual, 'Future'), get_time_period(ds_ssp585_seasonal, 'Future')
            ssp245_target, ssp245_target_seasonal = get_time_period(ds_ssp245_annual, 'Future'), get_time_period(ds_ssp245_seasonal, 'Future')
            G6sulfur_assmt, G6sulfur_assmt_seasonal = get_time_period(ds_G6sulfur_annual, 'Future'), get_time_period(ds_G6sulfur_seasonal, 'Future')
    
    
            
            process_and_save_maps(ssp245_baseline, ssp245_baseline_seasonal, 
                                  var=var, table='Amon', 
                                  label='SSP245_baseline', seasons=seasons)
            process_and_save_maps(ssp585_background, ssp585_background_seasonal, 
                                  var=var, table='Amon', 
                                  label='SSP585_background', seasons=seasons)
            process_and_save_maps(ssp245_target, ssp245_target_seasonal, 
                                  var=var, table='Amon', 
                                  label='SSP245_target', seasons=seasons)
            process_and_save_maps(G6sulfur_assmt, G6sulfur_assmt_seasonal, 
                                  var=var, table='Amon', 
                                  label='G6sulfur_assmt', seasons=seasons)
        except:
            print('{} failed, missing data?'.format(var))
        