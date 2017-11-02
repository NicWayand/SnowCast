
# coding: utf-8

# In[ ]:

get_ipython().magic(u'matplotlib inline')
#mpld3.enable_notebook()
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime
import xarray as xr
from astropy.io import ascii
import pytz
# OS interaction
import sys
import os
import glob
import wget
import seaborn as sns
sns.set_context("talk",font_scale=1.5)
sns.set_style('whitegrid')


# # User config

# In[ ]:




# In[ ]:

# Paths to user files
data_dir = os.path.normpath(r'F:\Work\e\Data\Obs\Canada_Project_Sites\CSAS_data') # Where to store data on local computer
git_dir  = os.path.normpath(r'C:\Users\new356\Google Drive\Python\SnowCast\In_Situ_Data') # This repo


# # Create paths
# 

# In[ ]:

# Data network
network = 'ABE_AGG_HIST'

# Location to download historical AB station data
download_dir = os.path.join(data_dir,network,'current')
# Make if does not exist
if not os.path.exists(download_dir):
    os.makedirs(download_dir)
    
# Netcdf file to save to
netcdf_dir   = os.path.join(data_dir,network,'netcdf')
# Make if does not exist
if not os.path.exists(netcdf_dir):
    os.makedirs(netcdf_dir)
netcdf_file_out =  os.path.join(netcdf_dir,'ABE_AGG_HIST.nc')

# # Metadata for AB pillows 
meta_file         = 'ABE_AGG_Station_Metadata.csv'
meta_file_path    = os.path.join(git_dir,'metadata',meta_file)


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:

os.chdir(download_dir)


# In[ ]:

sta_files = ['HourlyPrecipMatrix.csv',
            'HourlyRHMatrix.csv',
            'HourlyTaMatrix.csv',
            'HourlyUDMatrix.csv',
            'HourlyUSMatrix.csv']

Var_names = ['Precipitation','RealtiveHumidity','AirTemperature','WindDirection','WindSpeed']
Var_units = ['mm','%','C','m/s','degrees'] # not windspeed comes in as km/h, converted below


# In[ ]:

# Import each data variable
ds_dict   = {}
hrd_dict = {}
unit_dict = {}
for (i,cf) in enumerate(sta_files):
    print(cf)
    
    # Load in to python
    df_all = pd.read_csv(cf,index_col=0, engine='python', parse_dates=True, na_values=['---'])
    df_hdr_raw = df_all[0:14]
    df_dat = df_all[16:].convert_objects(convert_numeric=True) # Force to floats
    
    # Clean up messy header info
    df_hdr = df_hdr_raw.copy()
    for i2, row in df_hdr_raw.iterrows():
        df_hdr.loc[i2] = row.str.split(r'\t').str[1]
    
    # Fix dates
    df_dat.index = pd.DatetimeIndex(df_dat.index)
    
    # Drop rows with bad times (NAT)
    df_dat = df_dat[pd.notnull(df_dat.index)]
    
    var_name_full  = Var_names[i]
    var_units = Var_units[i]
    # Set column names as station ID
    df_dat.columns = df_hdr.loc['Station Number'].values
#     print(df_hdr.loc['Station Number'].values)
    # Set time
    df_dat.index.names = ['Time_MST'] # Kabir was unsure, check this time zone
    
    # Store as dict
    ds_dict[var_name_full]   = df_dat
    unit_dict[var_name_full] = var_units
    hrd_dict[var_name_full]  = df_hdr
    
# Memory clean up
df_all = None
df_hdr_raw = None
df_dat = None


# In[ ]:




# In[ ]:

# Constains duplicate time values...  So need to manualy create datasets, then fix time, then merge
da_fill_list = []
for k in ds_dict.keys():
    print(k)
    ds_temp = xr.DataArray(ds_dict[k], coords = {'Time_MST':ds_dict[k].index, 'staID':ds_dict[k].columns}, dims = ('Time_MST','staID'))
    ds_temp.name = k
    Time_new = np.arange(ds_temp.Time_MST.values[0], ds_temp.Time_MST.values[-1], dtype='datetime64[h]')
    ds_fill = ds_temp.reindex({'Time_MST':Time_new})
    da_fill_list.append(ds_fill)
ds_dict = None # Memory clean up


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:

# Merge into netcdf
ds = xr.merge(da_fill_list)
da_fill_list = None # Memory clean up


# In[ ]:

# Convert units
ds['WindSpeed'] = ds['WindSpeed'] * 1000/(60*60) # km/h to m/s


# In[ ]:

## ADD UNITS
# Add variable attributes (units), and fix variable names (remove spaces)
for cvar in ds.data_vars:
    # add units as attributes
    ds.get(cvar).attrs['unit'] = unit_dict[cvar]


# In[ ]:

# Grab metadata
metadata = hrd_dict['AirTemperature'] # Might need to merge all
metadata.columns = metadata.loc['Station Number'].values
metadata = metadata.transpose()


# In[ ]:




# In[ ]:

## Add station metadata
ds['station_name'] = xr.DataArray(metadata['Station Name'], coords={'staID':metadata.index}, dims='staID')
ds['Lat'] = xr.DataArray(metadata['Latitude'].astype(float),coords={'staID':metadata.index}, dims='staID')
ds['Lon'] = xr.DataArray(metadata['Longitude'].astype(float),coords={'staID':metadata.index}, dims='staID')
# Orig filese are missing Elevation, so used a separate metadata file
# Import metadata for each station
metadata_EL = pd.read_csv(meta_file_path,index_col='staID',delimiter=',',na_values=[0])
metadata_EL = metadata_EL.loc[ds.staID.values]
ds['Elevation'] = xr.DataArray(metadata_EL['Elevation'],coords={'staID':metadata_EL.index}, dims='staID')


# In[ ]:




# In[ ]:




# In[ ]:

# Move to coords
ds.set_coords(['station_name','Elevation','Lat','Lon'], inplace=True)
ds


# In[ ]:

# Find stations missing Elevation data (or is zero (wrong))
X = ds.Elevation.where(ds.Elevation.isnull(), drop=True)
list(zip(X.station_name.values,X.staID.values))


# In[ ]:

# # Reindex to have continous time steps
# Time_UTC_new = np.arange(ds.Time_UTC.values[0], ds.Time_UTC.values[-1], dtype='datetime64[h]')
# ds_fill = ds.reindex({'Time_UTC':Time_UTC_new})


# In[ ]:

# Add Network(s)
ds.coords['network'] = xr.DataArray([network for x in ds.staID], dims='staID')


# In[ ]:

# # Write out to netcdf by year
# os.chdir(netcdf_dir)
# years, datasets = zip(*ds.groupby('Time_MST.year'))
# paths = ['%s.nc' % y for y in years]
# xr.save_mfdataset(datasets, paths)


# In[ ]:

# Save as netcdf file
ds.to_netcdf(netcdf_file_out)


# In[ ]:



