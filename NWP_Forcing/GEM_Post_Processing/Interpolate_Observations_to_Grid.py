'''
Take netcdf file of observations of variable X, interpolate to a grid taken from a NWP model
'''
import xarray as xr
import pandas as pd
import numpy as np
import os
from scipy.interpolate import griddata

# Load in GEM model grid and point Observation data
ds_gem = xr.open_dataset(r'F:\Work\e\Data\Mesoscale\GEM\west\netcdf\GEM_2_5km_west_2017-03-23T01_00_00_000000000.nc')
print(ds_gem)
data_dir = r'F:\Work\e\Data\Obs\Canada_Project_Sites\CSAS_data'
obs_file = os.path.join(data_dir,'QC','Hourly_QC.nc')
ds_obs = xr.open_dataset(obs_file)
print(ds_obs)

# Define target grid
grid = 5
