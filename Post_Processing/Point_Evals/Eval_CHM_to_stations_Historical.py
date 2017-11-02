# Evaluate GEM-CHM Met and Snow over Historical Period

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
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import sys
import imp
import os
import seaborn as sns
plt.rcParams.update({'figure.max_open_warning': 0})

# General plotting settings
sns.set_style('ticks')
sns.set_context("talk", font_scale=3, rc={"lines.linewidth": 2.5})
fig_res = 90 # dpi

# Load in config file
#######  load user configurable paramters here    #######
# Check user defined configuraiton file
if len(sys.argv) != 3:
    sys.error('Requires two arguments [configuration file] [chm_run_dir]')

# Get name of configuration file/module
configfile = sys.argv[1]
chm_run_dir = str(sys.argv[2])

# Load in configuration file as module
X = imp.load_source('',configfile)

# Assign to local variables
data_dir = X.data_dir
git_dir   = X.git_dir

main_dir  = os.path.join(git_dir, 'CHM_Configs', chm_run_dir)
fig_dir   = os.path.join(main_dir , 'figures')

# Make fig dir
if not os.path.isdir(fig_dir):
    os.mkdir(fig_dir)

def save_figure(f,file_out,fig_res):
    f.savefig(file_out,bbox_inches='tight',dpi=fig_res)

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
#
# # Set font size
# font = {'weight' : 'bold',
#         'size'   : 24}
# matplotlib.rc('font', **font)

# Data files in
file_in = os.path.join(data_dir, 'QC', 'Hourly_QC.nc') # CRHO and other data
snow_survey_in  = os.path.join(data_dir, 'CRHO_HIST', 'netcdf', 'CRHO_Snow_Survey_Individual.nc')
EC_snow_course_in = os.path.join(data_dir, 'EC_Snow_Courses', 'netcdf', 'EC_Snow_Courses.nc')

# Load all obs
OBS_data = xr.open_dataset(file_in, engine='netcdf4') #.load()
# Rename obs variable names to model variable names
OBS_data.rename(vars_all, inplace=True);

# Filling in missing SW values at night (these were negative values that in QC SHOULD have been set to zero)
OBS_data['iswr'] = OBS_data['iswr'].fillna(0)
print('iswr fill is hack, need to fix upstream')

# Snow surveys
SS_data = xr.open_dataset(snow_survey_in,engine='netcdf4')
EC_data = xr.open_dataset(EC_snow_course_in)

# For current exp/folder, get netcdf file
c_mod_file = os.path.join(main_dir,'points','CHM_pts.nc')
Mod_data = xr.open_dataset(c_mod_file,engine='netcdf4')
dt_eval_hr = {'H':1}

EC_data.rename({'staID':'station', 'Time_UTC':'time', 'SnowDepth_point':'snowdepthavg', 'SWE_point':'swe'}, inplace=True);

# Find model cells that have forest on
# print(Mod_data.where(Mod_data.snow_load.sum(dim='time')>0, drop=True).station)
# print(Mod_data.where(Mod_data.snow_load.sum(dim='time')==0, drop=True).station.values)

# Function that makes two data sets common:
# Variables
# Time Step
# Stations
def make_common(ds_obs, ds_mod, dt_eval):
    # dt_eval = 'H' # MS (month start), H (hour), W (week)
    
    # Common variables
    obs_vars = ds_obs.data_vars
    mod_vars = ds_mod.data_vars
    # Find common variables
    com_vars = np.sort(list(set(obs_vars).intersection(mod_vars)))
    # Extract only common vars
    ds_obs=ds_obs[com_vars]
    ds_mod=ds_mod[com_vars]
    print("Common variables are:",com_vars)
    print("")

    # Common stations
    obs_sta = ds_obs['station'].values
    mod_sta = ds_mod['station'].values
    # Find common variables
    com_sta = np.sort(list(set(obs_sta).intersection(mod_sta)))
    print("Common stations are:",com_sta)
    print("")
    if not len(com_sta):
        print(obs_sta)
        print(mod_sta)
        raise ValueError('No common stations found')

    # Extract only common stations
    ds_obs=ds_obs.sel(station=com_sta)
    ds_mod=ds_mod.sel(station=com_sta)

    # Common time step
    # OBS
    obs_dt_in = (ds_obs.time.values[2]-ds_obs.time.values[1]).astype('timedelta64[h]').astype(float)
    if ((obs_dt_in != dt_eval_hr[dt_eval]) & (dt_eval!='exact')): # Check if we need to agg (exact skips agg (i.e. snow surveys))
        print('TODO: Current bug where if entire record is nans, it returns zeros...')
        print('Resampling OBS')
        obs_dt_val = ds_obs.resample(freq=dt_eval,dim='time',how='mean',label='left')
        if 'p' in ds_obs.data_vars:
            obs_dt_val['p'] = ds_obs['p'].resample(freq=dt_eval,dim='time',how='sum',label='left')

    else:
        obs_dt_val = ds_obs

    # Mod (advances by one hour if hourly in, hourly out... (why? orig label = left default??))
    mod_dt_in = (ds_mod.time.values[2]-ds_mod.time.values[1]).astype('timedelta64[h]').astype(float)
    if ((mod_dt_in != dt_eval_hr[dt_eval]) & (dt_eval!='exact')): # Check if we need to agg
        print('Resampling Model')
        mod_dt_val = ds_mod.resample(freq=dt_eval,dim='time',how='mean',label='left')
        if 'p' in ds_mod.data_vars:
            mod_dt_val['p'] = ds_mod['p'].resample(freq=dt_eval,dim='time',how='sum',label='left')

    else:
        mod_dt_val = ds_mod

    ## Common time period
    # agg time
    com_time = np.intersect1d(mod_dt_val.time,obs_dt_val.time)
    print("Common time is:",com_time[0], " to ", com_time[-1])
    print("")
    obs_dt_val = obs_dt_val.sel(time=com_time).load()
    mod_dt_val = mod_dt_val.sel(time=com_time).load()

#     # Clean up
#     ds_mod = None
#     ds_obs = None
    
    return (obs_dt_val, mod_dt_val)

# Get common obs and model
(obs_dt_val, mod_dt_val) = make_common(OBS_data, Mod_data, 'H')

# Get common survey and model
# (obs_Survey, mod_Survey) = make_common(EC_data, Mod_data, 'exact')

# Memory Clean up
OBS_data = None
Mod_data = None


# ## Calculate Stats

# def calc_stats(ds_obs, ds_mod):
#     # RMSE
#     rmse =

# Plot settings
import seaborn as sns
sns.set_style('whitegrid')
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})

sta_list = np.sort(mod_dt_val.station)

## Time series plots
#cmap = colors.ListedColormap(['white', 'red'])
sta_clrs = ['b','g','r','c','m','y','k', 'b','g','r','c','m','y','k', 'b','g','r','c','k', 'b','g','r',
            'c','m','y','k', 'b','g','r','c','k', 'b','g','r','c','m','y','k', 'b','g','r','c','k']
sta_cmap = dict(zip(mod_dt_val.station.values, sta_clrs))

obs_linewidth = 2.5
mod_linewidth = 2


# How well can we Simulate point SWE?
Vars_to_plot = ['swe']
# Create a new figure and subplots for each variable
(f, ax1) = plt.subplots(1, 1, sharey=False)
f.set_size_inches(16, 6)
#### Loop through each station
v_c = 0
for cvar in Vars_to_plot:
    #     h_m = [] # init handles for plots
    #     h_o = []
    #### Loop through each station
    
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
file_out = os.path.join(fig_dir, 'Point_' + Vars_to_plot[0] + '.png')
save_figure(f,file_out,fig_res)


#
# for csta in mod_dt_val.station:
#     plt.figure()
#     plt.title(str(csta.values))
#     plt.plot(mod_dt_val.time, mod_dt_val['swe'].sel(station=csta))

Vars_to_plot = ['snowdepthavg']

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
file_out = os.path.join(fig_dir, 'Point_' + Vars_to_plot[0] + '.png')
save_figure(f,file_out,fig_res)



Vars_to_plot = ['t']
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
file_out = os.path.join(fig_dir, 'Point_' + Vars_to_plot[0] + '.png')
save_figure(f,file_out,fig_res)

# plt.show()