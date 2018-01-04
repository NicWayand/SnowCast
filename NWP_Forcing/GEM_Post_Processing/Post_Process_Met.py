'''
Take netcdf file of observations of variable X, interpolate to a grid taken from a NWP model
'''
import xarray as xr
import pandas as pd
import numpy as np
import os
import sys
import datetime
from scipy.interpolate import griddata
import imp
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import PP
import idw
import json

# General plotting settings
sns.set_style('whitegrid')
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
fig_res = 90 # dpi

# Load in config file
#######  load user configurable paramters here    #######
# Check user defined configuraiton file
if len(sys.argv) != 2:
    sys.exit('Requires 1 argument [configuration file] ')

# Get name of configuration file/module
configfile = sys.argv[1]
# Load in configuration file as module
X = imp.load_source('',configfile)

# Assign to local variables
data_dir = X.data_dir
git_dir = X.git_dir

# Interpolated obs file out
obs_file = os.path.join(data_dir, 'QC', 'Gridded.nc')
ds_obs = xr.open_dataset(obs_file)

# Load in info on GEM data (in ascii format)
main_gem_dir = '/media/data3/nicway/GEM/archive/SOAP/'
ascii_in_dir = 'ascii_HRDPS_SnowCast_full'
GEM_config_file = os.path.join(main_gem_dir, ascii_in_dir, 'SOAP_forcing.json')
gem_files = json.load(open(GEM_config_file))

# Output corrected GEM dir
ascii_out_dir = 'precip_test'
output_gem_dir = os.path.join(main_gem_dir, ascii_out_dir)
GEM_config_file_out = os.path.join(main_gem_dir, ascii_out_dir, 'SOAP_forcing.json')
if not os.path.exists(output_gem_dir):
    os.makedirs(output_gem_dir)



# Ensemble number
# N_e = 8

def naive_fast(latvar,lonvar,lat0,lon0):
    # Read latitude and longitude from file into numpy arrays
    latvals = latvar[:]
    lonvals = lonvar[:]
    dist_sq = (latvals-lat0)**2 + (lonvals-lon0)**2
    minindex_flattened = dist_sq.argmin()  # 1D index of min element
    iy_min,ix_min = np.unravel_index(minindex_flattened, latvals.shape)
    return iy_min,ix_min

# Initialize df to save stats
df_stats = pd.DataFrame(index=gem_files.keys(),
                        columns=['orig_Bias', 'adj_Bias'])

# Loop through each GEM file (grid cell)
for cpt, cinfo in gem_files.iteritems():
    print cpt
    cFile = open(cinfo['file'],'r')
    cLat = cinfo['latitude']
    cLon = cinfo['longitude']
    df = pd.read_csv(cFile, sep="\t", parse_dates=True)
    cFile.close()
    df.set_index('datetime', inplace=True)
    df.index = pd.to_datetime(df.index)
    p_mod = df['p'] # get series
    p_mod.name = 'model'

    # Find nearest obs cell
    (c_y, c_x) = naive_fast(ds_obs.gridlat_0.values,
                            ds_obs.gridlon_0.values,
                            cLat, cLon)
    p_obs = ds_obs.isel(xgrid_0=c_x, ygrid_0=c_y).IncrementalPrecipitationA.to_series()
    p_obs = p_obs*1000 # m to mm
    p_obs.name = 'observed'

    # Trim to common period
    df_mrg = pd.concat([p_mod, p_obs], axis=1, join='inner')

    # Test plots
    # df_mrg.plot()
    # df_mrg.cumsum().plot()
    # plt.figure()
    # plt.scatter(df_mrg['observed'],df_mrg['model'])

    # Select training period (year/season/month)
    t_start = '2014-10-01'
    t_end = '2015-09-30'
    df_t = df_mrg.loc[t_start:t_end]

    # Fit model
    # https://stackoverflow.com/questions/39330483/fitting-location-parameter-in-the-gamma-distribution-with-scipy
    t_sum = df_t.dropna().sum(axis=0) # Accumluated precip in t period (mm)
    cBiasRatio = t_sum['observed']/t_sum['model']
    gem_files[cpt]['pBiasRatio'] = cBiasRatio # store in json dictionary

    # Evaluate PP model to obs over evaluation period
    e_start = '2015-10-01'
    e_end = '2016-09-30'
    df_eval = df_mrg.loc[e_start:e_end]
    # Apply model to evaluation period
    df_eval['model_adj'] = df_eval.model * cBiasRatio

    # Evaluate fit
    e_sum = df_eval.dropna().sum(axis=0) # Accumluated precip in t period (mm)
    new_Bias = e_sum['observed']/e_sum['model_adj']
    orig_Bias = e_sum['observed']/e_sum['model']
    df_stats.set_value(cpt, 'orig_Bias', orig_Bias)
    df_stats.set_value(cpt, 'adj_Bias', new_Bias)

    # Apply to full GEM/model period and export to new ascii forcing dir
    df['p'] = df.p * cBiasRatio
    file_out = open(os.path.join(output_gem_dir, cpt+'.chm'),'w')
    df.to_csv(file_out,sep='\t',date_format='%Y%m%dT%H%M%S')
    file_out.close()
    # Update config file
    gem_files[cpt]['file'] = file_out

# Write out json file with metadata for these points
with open(GEM_config_file_out, 'w') as outfile:
    json.dump(gem_files, outfile,indent=4, sort_keys=True)

# Summary results of all obs grid points
df_stats.plot()
plt.show()

# x = 5

# Apply

# Export PP model to original ascii format
