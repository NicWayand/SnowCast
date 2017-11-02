
# coding: utf-8

# # SnowCast Spatial Error Maps

# In[10]:

# Import modules
# ipython magic to plot in line
get_ipython().magic('matplotlib inline')
#import mpld3
#mpld3.enable_notebook()
import matplotlib
#matplotlib.style.use('ggplot')
import numpy as np
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from astropy.io import ascii
import sklearn.metrics as skm
import pytz
from sklearn.metrics import mean_squared_error
from math import sqrt
import scipy
# OS interaction
import sys
import os
import glob
import time
plt.rcParams.update({'figure.max_open_warning': 0})
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import xarray as xr
import seaborn as sns


# In[ ]:




# In[11]:

# Config info
#Stations_all = ['BNS','CRG','CRN','FLG','FRG','FRS','PWL']
#Station_full_names = ['Bonsai','Can. Ridge','Can. Ridge North','Fortress Ledge','Fortress Ridge','Fortress Ridge South','Powerline']
# CRHO
#Stations_all = ['BNS', 'BRP', 'BWH', 'CNT', 'CRG', 'CRN', 'FLG', 'FRG', 'FRS',
#       'FSR', 'HLN', 'HMW', 'LLF', 'PWL', 'PYT', 'SIB', 'UPC', 'UPF', 'VVW']

# Make dictionary of obs:model variable names to compare
# model:obs
vars_all = {'AirtemperatureA':'t','AirMoistureContentA':'rh','IncrementalPrecipitationA':'p',
            'ScalarWindSpeedA':'U_2m_above_srf','DownwardSolarRadiation':'iswr','DownwardTerrestrialRad':'ilwr',
            'UpwardTerrestrialRad':'ilwr_out',
            'SnowWaterEquivelentA':'swe','SnowDepthA':'snowdepthavg','WindDirectionatA':'vw_dir',
           'Time_UTC':'time','staID':'station'}

plot_key = {'ilwr_out':'Outgoing Longwave','T_s_0':'Surface temperature','t':'Air Temperature','rh':'Relative Humidity',
            'p':'Precipitation','ilwr':'Incoming Longwave','iswr':'Shortwave Radiation',
            'U_2m_above_srf':'Wind Speed','vw_dir':'Wind Direction','swe':'Snow Water Equivalent','snowdepthavg':'Snowdepth'}

ylabel_unit = {'ilwr_out':'W m-2','G':'W m-2','T_s_0':'(C)','t':'C','rh':'%','p':'m','ilwr':'W m-2','iswr':'W m-2',
            'U_2m_above_srf':'m/s','vw_dir':'degrees true north','swe':'m','snowdepthavg':'m'}




# DS = datetime.datetime(2014,11,13) #2014-11-12T17
# DE = datetime.datetime(2016,5,1)

cExp=0


# In[12]:

EXP_names = ['HRDPS_Historical']
# EXP_names = ['forecast_CRHO_spinup','GDPS_Current']


# In[13]:

# Set font size
font = {'weight' : 'bold',
        'size'   : 24}
matplotlib.rc('font', **font)


# In[14]:

# CHM output
mod_dir   = os.path.normpath(r'C:\\Users\new356\Model_Output\CHM\SnowCast')

nc_file_out = 'CHM_pts.nc'
dt_eval = 'MS' # MS (month start), H (hour), W (week)
dt_eval_hr = {'H':1,'MS':999999}

# For current exp/folder, get netcdf file
c_mod_file = os.path.join(mod_dir, EXP_names[cExp],'points',nc_file_out)
Mod_data = xr.open_dataset(c_mod_file,engine='netcdf4')


# In[ ]:




# In[ ]:




# In[ ]:




# In[15]:

# Load in Provences shapefile
p_file = r'F:\Work\e\Data\Shape_files\CAN_adm1.shp'
p_sh = list(shpreader.Reader(p_file).geometries())


# In[16]:

## Observed
file_in = os.path.normpath(r'F:\Work\e\Data\Obs\Canada_Project_Sites\CSAS_data\QC\Hourly_QC.nc') # CRHO and other data
# Load all obs
CRHO_data = xr.open_dataset(file_in,engine='netcdf4') #.load()
# Rename obs variable names to model variable names
CRHO_data.rename(vars_all,inplace=True);
#CRHO_data['G'].sel(station='CRN') = CRHO_data['G'].sel(station='CRN') * -1
# Filling in missing SW values at night (these were negative values that in QC SHOULD have been set to zero)
CRHO_data['iswr'] = CRHO_data['iswr'].fillna(0)


# In[ ]:




# In[17]:

# Load in dem
dem_file = r'F:\Work\e\Data\DEMs\NA.tif'
dem = xr.open_rasterio(dem_file).sel(band=1).drop('band')



# In[ ]:




# In[ ]:




# In[20]:

cmap_elev = matlibplot.colors.ListedColormap(sns.color_palette("Greys", 10))
cmap_network = {'bcRiverForecastCenter':'r', 'environmentAlberta':'k','CRHO':'b','ABE_AGG_HIST':'g'}


# In[ ]:

sns.set_style('whitegrid')
sns.set_context("talk", font_scale=2, rc={"lines.linewidth": 2})


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:

# t = CRHO_data.t.isel(station=15).sel(time=slice('2014-08-01','2014-08-02'))
# t = t.where(t.notnull(),drop=True)
# t.plot()


# In[ ]:




# In[ ]:

# Find model cells that have forest on
Mod_data.where(Mod_data.snow_load.sum(dim='time')>0, drop=True).station


# In[ ]:

Mod_data.where(Mod_data.snow_load.sum(dim='time')==0, drop=True).station.values


# In[ ]:




# In[ ]:

# Common variables
obs_vars = CRHO_data.data_vars
mod_vars = Mod_data.data_vars
# Find common variables
com_vars = np.sort(list(set(obs_vars).intersection(mod_vars)))
# Extract only common vars
CRHO_data=CRHO_data[com_vars]
Mod_data=Mod_data[com_vars]

# Common stations
obs_sta = CRHO_data['station'].values
mod_sta = Mod_data['station'].values
# Find common variables
com_sta = np.sort(list(set(obs_sta).intersection(mod_sta)))
print(com_sta)
# Extract only common stations
CRHO_data=CRHO_data.sel(station=com_sta)
Mod_data=Mod_data.sel(station=com_sta)

# TODO: Should we force model to be missing where obs are missing? To make weights of means equal? 

# Common time step
# OBS
obs_dt_in = (CRHO_data.time.values[2]-CRHO_data.time.values[1]).astype('timedelta64[h]').astype(float)
if (obs_dt_in != dt_eval_hr[dt_eval]): # Check if we need to agg
    print('TODO: Current bug where if entire record is nans, it returns zeros...')
    print('Resampling OBS')
    obs_dt_val = CRHO_data.resample(freq=dt_eval,dim='time',how='mean',label='left')
    if 'p' in CRHO_data.data_vars:
        obs_dt_val['p'] = CRHO_data['p'].resample(freq=dt_eval,dim='time',how='sum',label='left')
    
else:
    obs_dt_val = CRHO_data

# Mod (advances by one hour if hourly in, hourly out... (why? orig label = left default??))
mod_dt_in = (Mod_data.time.values[2]-Mod_data.time.values[1]).astype('timedelta64[h]').astype(float)
if (mod_dt_in != dt_eval_hr[dt_eval]): # Check if we need to agg
    print('Resampling Model')
    mod_dt_val = Mod_data.resample(freq=dt_eval,dim='time',how='mean',label='left')
    if 'p' in Mod_data.data_vars:
        mod_dt_val['p'] = Mod_data['p'].resample(freq=dt_eval,dim='time',how='sum',label='left')
    
else:
    mod_dt_val = Mod_data

## Common time period
# agg time
com_time = np.intersect1d(mod_dt_val.time,obs_dt_val.time)
print(com_time[0], " to ", com_time[-1])
obs_dt_val = obs_dt_val.sel(time=com_time).load()
mod_dt_val = mod_dt_val.sel(time=com_time).load()

# Clean up
Mod_data = None
CRHO_data = None


# In[19]:

obs_dt_val


# In[ ]:

lat_r = ds.Lat.max()-ds.Lat.min()
lon_r = ds.Lon.max()-ds.Lon.min()
bdy = 0.2
box = [ds.Lon.min()-lon_r*bdy, ds.Lon.max()+lon_r*bdy, ds.Lat.min()-lat_r*bdy, ds.Lat.max()+lat_r*bdy]
# box = [ds.Lon.min(), ds.Lon.max(), ds.Lat.min(), ds.Lat.max()]


# In[ ]:

# Set nan in dem
dem = dem.where(dem>0)
dem = dem.where((dem.x>box[0]) & (dem.x<box[1]) & (dem.y>box[2]) & (dem.y<box[3]), drop=True)


# In[ ]:

fig = plt.figure(figsize=(20, 20))
ax1 = plt.axes(projection=ccrs.AlbersEqualArea())
# ax1.set_extent(box)
ax1.imshow(np.flipud(dem.values), extent=[np.min(dem.x), np.max(dem.x),
                                             np.min(dem.y), np.max(dem.y)], aspect=ax1.get_aspect())
# ax1.set_title('Elevation')
for c_net in set(ds.network.values):
    lat_pts = ds.Lat.sel(staID=(ds.where(ds.network==c_net, drop=True).network).staID).values
    lon_pts = ds.Lon.sel(staID=(ds.where(ds.network==c_net, drop=True).network).staID).values
    I_not_nan = ~np.isnan(lat_pts) & ~np.isnan(lon_pts)
    lat_pts = lat_pts[I_not_nan]
    lon_pts = lon_pts[I_not_nan]
    
    ax1.scatter(lon_pts, lat_pts, transform=ccrs.AlbersEqualArea(), s=50, c=cmap_network[c_net], zorder=100, label=c_net) #yc, xc -- lists or numpy arrays

ax1.add_geometries(p_sh, ccrs.AlbersEqualArea(),
                  edgecolor='black', facecolor='none', alpha=0.5)
plt.legend()


# In[ ]:




# In[ ]:




# In[ ]:



