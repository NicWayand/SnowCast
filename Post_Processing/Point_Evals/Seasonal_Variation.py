
# coding: utf-8

# # Evaluate GEM-CHM Met and Snow over Historical Period

# ## Point stations come from
# ## 1) CRHO (many vars)
# ## 2) BC Hydro (swe, SD, T, P)
# ## 3) AB snow pillows (swe, SD)

# ## Time Periods (UTC)
#     ## CHM:   Nov 2014 to Aug 2017
#     ## CRHO:      2002 to May 2016
#     ## BC:        2016 to 2017 WY (Missing historical data)   
#     ## AB:        1984 to Present
# ## Common time period
#     ## Nov 2014 to Aug 2017

# In[ ]:

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


# In[ ]:




# In[ ]:

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


# In[ ]:

EXP_names = ['HRDPS_Historical']
# EXP_names = ['forecast_CRHO_spinup','GDPS_Current']


# In[ ]:

# Set font size
font = {'weight' : 'bold',
        'size'   : 24}
matplotlib.rc('font', **font)


# In[ ]:

# CHM output
mod_dir   = os.path.normpath(r'C:\\Users\new356\Model_Output\CHM\SnowCast')

nc_file_out = 'CHM_pts.nc'
dt_eval = 'MS' # MS (month start), H (hour), W (week)
dt_eval_hr = {'H':1,'MS':999999}
## Observed
file_in = os.path.normpath(r'F:\Work\e\Data\Obs\Canada_Project_Sites\CSAS_data\QC\Hourly_QC.nc') # CRHO and other data


# In[ ]:




# In[ ]:

# Load all obs
CRHO_data = xr.open_dataset(file_in,engine='netcdf4') #.load()
# Rename obs variable names to model variable names
CRHO_data.rename(vars_all,inplace=True);
# Flip direction of ground heat flux to match model (postive upward into snowpack)
# CRHO_data['G'] = CRHO_data['G'] * -1
# CRN ground heat flux wrong sign
#CRHO_data['G'].sel(station='CRN') = CRHO_data['G'].sel(station='CRN') * -1
# Filling in missing SW values at night (these were negative values that in QC SHOULD have been set to zero)
CRHO_data['iswr'] = CRHO_data['iswr'].fillna(0)


# In[ ]:

# For current exp/folder, get netcdf file
c_mod_file = os.path.join(mod_dir, EXP_names[cExp],'points',nc_file_out)
Mod_data = xr.open_dataset(c_mod_file,engine='netcdf4')


# In[ ]:




# In[ ]:

# for csta in Mod_data.station:
#     plt.figure()
#     plt.title(str(csta.values))
#     plt.plot(Mod_data.time, Mod_data['iswr_subcanopy'].sel(station=csta))


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


# In[ ]:




# In[ ]:




# ## Calculate Stats

# In[ ]:

# def calc_stats(ds_obs, ds_mod):
#     # RMSE
#     rmse = 


# In[ ]:

# PLot settings
import seaborn as sns
sns.set_style('whitegrid')
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})


# In[ ]:




# In[ ]:




# ## How well can we Simulate point SWE?

# In[ ]:

Vars_to_plot = ['swe']

## Time series plots
#cmap = colors.ListedColormap(['white', 'red'])
sta_clrs = ['b','g','r','c','m','y','k', 'b','g','r','c','m','y','k', 'b','g','r','c','k', 'b','g','r',
            'c','m','y','k', 'b','g','r','c','k', 'b','g','r','c','m','y','k', 'b','g','r','c','k']
sta_cmap = dict(zip(com_sta, sta_clrs))
#date_range = [datetime.datetime(2014,11,13),datetime.datetime(2014,12,4)]

obs_linewidth = 2.5
mod_linewidth = 2

# Create a new figure and subplots for each variable
(f, ax1) = plt.subplots(1, 1, sharey=False)
f.set_size_inches(16, 6)
# ax1 = ax1.flatten()

#### Loop through each station
v_c = 0
for cvar in Vars_to_plot:
    #     h_m = [] # init handles for plots
    #     h_o = []
    #### Loop through each station
    sta_list = np.sort(com_sta)
    for csta in sta_list:
        # Get mod and obs arrays
        c_mod = mod_dt_val[cvar].sel(station=csta)
        c_obs = obs_dt_val[cvar].sel(station=csta)
        # Find non-nan values
        I_not_nans = ~c_mod.isnull() & ~c_obs.isnull()
        # If any non-nan, take metrics
        if np.any(I_not_nans)  and sum(c_obs)!=0: # last check is to remove SW stations with all zeros 
            # print ('CHECK THIS is not shifting the time...')
            print(csta)
            c_mod_g=c_mod #[I_not_nans]
            c_obs_g=c_obs #[I_not_nans]
            
            # Plot model/obs to get quick check
            if(cvar!='p'):
                c_mod_g.plot.line(color=sta_cmap[csta],ax=ax1,linestyle='--',linewidth=mod_linewidth,label='_nolegend_') #marker='.',
                c_obs_g.plot.line(color=sta_cmap[csta],linewidth=obs_linewidth,ax=ax1,label=str(obs_dt_val.station_name.sel(station=csta).values))
            else:
                # Check there is any precip measured (bug from above where missing all nanas, are avveraged to zeros...)
                if np.all(sum(c_mod_g)>0 and sum(c_obs_g)>0):
                    np.cumsum(c_mod_g).plot.line(color=sta_cmap[csta],ax=ax1,linestyle='--',linewidth=mod_linewidth) #marker='.',
                    np.cumsum(c_obs_g).plot.line(color=sta_cmap[csta],linewidth=obs_linewidth,ax=ax1)
            ax1.set_title(plot_key[cvar])
            ax1.set_ylabel('SWE ('+ylabel_unit[cvar]+')')
            
            
    # incremetn axes                        
    v_c = v_c + 1
f.tight_layout()
plt.legend(loc='upper left')


# In[ ]:

# mod_dt_val['swe'].sel(station='PWL').values[0:1000]


# In[ ]:




# In[ ]:

Vars_to_plot = ['snowdepthavg']

## Time series plots
#cmap = colors.ListedColormap(['white', 'red'])
sta_clrs = ['b','g','r','c','m','y','k', 'b','g','r','c','m','y','k', 'b','g','r','c','m']
sta_cmap = dict(zip(com_sta, sta_clrs))
#date_range = [datetime.datetime(2014,11,13),datetime.datetime(2014,12,4)]

obs_linewidth = 2
mod_linewidth = 2

# Create a new figure and subplots for each variable
(f, ax1) = plt.subplots(1, 1, sharey=False)
f.set_size_inches(16, 6)
# ax1 = ax1.flatten()

#### Loop through each station
v_c = 0
for cvar in Vars_to_plot:
    #     h_m = [] # init handles for plots
    #     h_o = []
    #### Loop through each station
    sta_list = np.sort(com_sta)
    for csta in sta_list:
        # Get mod and obs arrays
        c_mod = mod_dt_val[cvar].sel(station=csta)
        c_obs = obs_dt_val[cvar].sel(station=csta)
        # Find non-nan values
        I_not_nans = ~c_mod.isnull() & ~c_obs.isnull()
        # If any non-nan, take metrics
        if np.any(I_not_nans)  and sum(c_obs)!=0: # last check is to remove SW stations with all zeros 
            # print ('CHECK THIS is not shifting the time...')
            print(csta, " ", str(obs_dt_val.station_name.sel(station=csta).values))
            c_mod_g=c_mod #[I_not_nans]
            c_obs_g=c_obs #[I_not_nans]
            
            # Plot model/obs to get quick check
            if(cvar!='p'):
                c_mod_g.plot.line(marker='.',color='r',ax=ax1,linestyle='--',linewidth=mod_linewidth)
                c_obs_g.plot.line(color='b',linewidth=obs_linewidth,ax=ax1)
            else:
                # Check there is any precip measured (bug from above where missing all nanas, are avveraged to zeros...)
                if np.all(sum(c_mod_g)>0 and sum(c_obs_g)>0):
                    np.cumsum(c_mod_g).plot.line(marker='.',color='r',ax=ax1,linestyle='--',linewidth=mod_linewidth)
                    np.cumsum(c_obs_g).plot.line(color='b',linewidth=obs_linewidth,ax=ax1)
            ax1.set_title(plot_key[cvar])
            ax1.set_ylabel('Snow Depth ('+ylabel_unit[cvar]+')')
            
            
    # incremetn axes                        
    v_c = v_c + 1
f.tight_layout()


# In[ ]:

Vars_to_plot = ['t']

## Time series plots
#cmap = colors.ListedColormap(['white', 'red'])
sta_clrs = ['b','g','r','c','m','y','k', 'b','g','r','c','m','y','k', 'b','g','r','c','m']
sta_cmap = dict(zip(com_sta, sta_clrs))
#date_range = [datetime.datetime(2014,11,13),datetime.datetime(2014,12,4)]

obs_linewidth = 2
mod_linewidth = 2

# Create a new figure and subplots for each variable
(f, ax1) = plt.subplots(1, 1, sharey=False)
f.set_size_inches(16, 6)
# ax1 = ax1.flatten()

#### Loop through each station
v_c = 0
for cvar in Vars_to_plot:
    #     h_m = [] # init handles for plots
    #     h_o = []
    #### Loop through each station
    sta_list = np.sort(com_sta)
    for csta in sta_list:
        # Get mod and obs arrays
        c_mod = mod_dt_val[cvar].sel(station=csta)
        c_obs = obs_dt_val[cvar].sel(station=csta)
        # Find non-nan values
        I_not_nans = ~c_mod.isnull() & ~c_obs.isnull()
        # If any non-nan, take metrics
        if np.any(I_not_nans)  and sum(c_obs)!=0: # last check is to remove SW stations with all zeros 
            # print ('CHECK THIS is not shifting the time...')
            print(csta)
            c_mod_g=c_mod #[I_not_nans]
            c_obs_g=c_obs #[I_not_nans]
            
            # Plot model/obs to get quick check
            if(cvar!='p'):
                c_mod_g.plot.line(marker='.',color='r',ax=ax1,linestyle='--',linewidth=mod_linewidth)
                c_obs_g.plot.line(color='b',linewidth=obs_linewidth,ax=ax1)
            else:
                # Check there is any precip measured (bug from above where missing all nanas, are avveraged to zeros...)
                if np.all(sum(c_mod_g)>0 and sum(c_obs_g)>0):
                    np.cumsum(c_mod_g).plot.line(marker='.',color='r',ax=ax1,linestyle='--',linewidth=mod_linewidth)
                    np.cumsum(c_obs_g).plot.line(color='b',linewidth=obs_linewidth,ax=ax1)
            ax1.set_title(plot_key[cvar])
            ax1.set_ylabel('SWE ('+ylabel_unit[cvar]+')')
            
            
    # incremetn axes                        
    v_c = v_c + 1
f.tight_layout()


# In[ ]:

# How much of monthly mean do we allow to be missing? (Need to screen obs and model so they have the same missing periodos)


# In[ ]:

Vars_to_plot = ['t','rh','U_2m_above_srf','p','ilwr','iswr']

## Time series plots
#cmap = colors.ListedColormap(['white', 'red'])
sta_clrs = ['b','g','r','c','m','y','k', 'b','g','r','c','m','y','k', 'b','g','r','c','m']
sta_cmap = dict(zip(com_sta, sta_clrs))
#date_range = [datetime.datetime(2014,11,13),datetime.datetime(2014,12,4)]

obs_linewidth = 2
mod_linewidth = 2

# Create a new figure and subplots for each variable
(f, ax1) = plt.subplots(2, 3, sharey=False)
f.set_size_inches(16, 8)
ax1 = ax1.flatten()

#### Loop through each station
v_c = 0
for cvar in Vars_to_plot:
    #     h_m = [] # init handles for plots
    #     h_o = []
    #### Loop through each station
    sta_list = np.sort(com_sta)
    for csta in sta_list:
        # Get mod and obs arrays
        c_mod = mod_dt_val[cvar].sel(station=csta)
        c_obs = obs_dt_val[cvar].sel(station=csta)
        #Find non-nan values
        I_not_nans = ~c_mod.isnull() & ~c_obs.isnull()
        # If any non-nan, take metrics
        if np.any(I_not_nans)  and sum(c_obs)!=0: # last check is to remove SW stations with all zeros 
            # print ('CHECK THIS is not shifting the time...')
            c_mod_g=c_mod #[I_not_nans]
            c_obs_g=c_obs #[I_not_nans]
            
            # Plot model/obs to get quick check
            if(cvar!='p'):
                c_mod_g.plot.line(color='r',ax=ax1[v_c],linestyle='--',linewidth=mod_linewidth)
                c_obs_g.plot.line(color='b',linewidth=obs_linewidth,ax=ax1[v_c])
            else:
                # Check there is any precip measured (bug from above where missing all nanas, are avveraged to zeros...)
                if np.all(sum(c_mod_g)>0 and sum(c_obs_g)>0):
                    c_mod_g.cumsum(dim='time').plot.line(color='r',ax=ax1[v_c],linestyle='--',linewidth=mod_linewidth)
                    c_obs_g.cumsum(dim='time').plot.line(color='b',linewidth=obs_linewidth,ax=ax1[v_c])
            ax1[v_c].set_title(plot_key[cvar])
            ax1[v_c].set_ylabel(ylabel_unit[cvar])
            
            
    # incremetn axes                        
    v_c = v_c + 1
f.tight_layout()


# In[ ]:




# In[ ]:




# In[ ]:

# mod_clean = c_mod_g.where(I_not_nans, drop=True).values
# obs_clean = c_obs_g.where(I_not_nans, drop=True).values

# print(cvar)
# print(np.mean(mod_clean)-np.mean(obs_clean))
# print(sqrt(mean_squared_error(mod_clean, obs_clean)))


# In[ ]:




# In[ ]:

# plt.figure()
# plt.plot(c_mod,c_obs,'ko')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.ylim([0,1000])
# plt.xlim([0,1000])


# In[ ]:

# import mpld3
# mpld3.enable_notebook()
# plt.figure()
# c_mod_g.plot.line(color='r',ax=plt.gca())
# c_obs_g.plot.line(color='b',ax=plt.gca())


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



