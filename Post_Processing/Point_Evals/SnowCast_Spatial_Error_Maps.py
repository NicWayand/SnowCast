# Standard modules
import matplotlib as mpl
# mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import sys
import imp
import os
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import seaborn as sns
plt.rcParams.update({'figure.max_open_warning': 0})
# SnowCast modules
import CHM_functions as chmF

# General plotting settings
sns.set_style('whitegrid')
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})
fig_res = 90 # dpi

# Load in config file
#######  load user configurable paramters here    #######
# Check user defined configuraiton file
if len(sys.argv) != 3:
    sys.exit('Requires two arguments [configuration file] [chm_run_dir]')

# Get name of configuration file/module
configfile = sys.argv[1]
chm_run_dir = str(sys.argv[2])

if chm_run_dir=='forecast_CRHO_spinup':
    c_run_dt_in = 'H'
elif chm_run_dir=='HRDPS_Current_BS':
    c_run_dt_in = 'H'
elif chm_run_dir=='HRDPS_Historical':
    c_run_dt_in = 'H'
elif chm_run_dir=='GDPS_Current':
    c_run_dt_in = '3H'
else:
    sys.exit('Model run name not found')

# Load in configuration file as module
X = imp.load_source('',configfile)

# Assign to local variables
data_dir = X.data_dir
git_dir = X.git_dir

main_dir = os.path.join(git_dir, 'CHM_Configs', chm_run_dir)
fig_dir = os.path.join(main_dir , 'figures', 'Error_Maps')

# Make fig dir
if not os.path.isdir(os.path.join(main_dir, 'figures')):
    os.mkdir(os.path.join(main_dir, 'figures'))
if not os.path.isdir(fig_dir):
    os.mkdir(fig_dir)



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

ylabel_unit = {'ilwr_out':'W m-2','G':'W m-2','T_s_0':'°C','t':'°C','rh':'%','p':'m','ilwr':'W m-2','iswr':'W m-2',
            'U_2m_above_srf':'m/s','vw_dir':'degrees true north','swe':'m','snowdepthavg':'m'}

###################################
# Load data in
###################################

# File paths
file_in = os.path.join(data_dir, 'QC', 'Hourly_QC.nc') # CRHO and other data
snow_survey_in  = os.path.join(data_dir, 'CRHO_HIST', 'netcdf', 'CRHO_Snow_Survey_Individual.nc')
EC_snow_course_in = os.path.join(data_dir, 'EC_Snow_Courses', 'netcdf', 'EC_Snow_Courses.nc')
dem_file = os.path.join(data_dir, 'Static_data', 'SnowCast.tif')
p_file = os.path.join(data_dir, 'Static_data', 'CAN_adm1.shp')
c_mod_file = os.path.join(main_dir,'points','CHM_pts.nc')

# Load all obs
OBS_data = xr.open_dataset(file_in, engine='netcdf4') #.load()

# Snow surveys
SS_data = xr.open_dataset(snow_survey_in,engine='netcdf4')
EC_data = xr.open_dataset(EC_snow_course_in)

# For current exp/folder, get netcdf file
Mod_data = xr.open_dataset(c_mod_file,engine='netcdf4')

# Load in dem
dem = xr.open_rasterio(dem_file).sel(band=1).drop('band')
# Provences
p_sh = list(shpreader.Reader(p_file).geometries())

####################################
# Modify data
####################################
# Rename obs variable names to model variable names
OBS_data.rename(vars_all, inplace=True);
EC_data.rename({'staID':'station', 'Time_UTC':'time', 'SnowDepth_point':'snowdepthavg', 'SWE_point':'swe'}, inplace=True);
# Filling in missing SW values at night (these were negative values that in QC SHOULD have been set to zero)
# OBS_data['iswr'] = OBS_data['iswr'].fillna(0)
# print('iswr fill is hack, need to fix upstream')

dt_eval_hr = {'H':1, '3H':3, 'MS':999999, 'W':999999} # This converts resample() strs to int hours. Use 999 if N/A.

#################################################
# Resample and make common
#################################################

# Get common obs and model
(obs_dt_val, mod_dt_val) = chmF.make_common(OBS_data, Mod_data, c_run_dt_in, dt_eval_hr)

# Memory Clean up
OBS_data = None
Mod_data = None


#############################################################
# Set up plotting info
#############################################################
sta_list = np.sort(mod_dt_val.station)

lat_r = obs_dt_val.Lat.max()-obs_dt_val.Lat.min()
lon_r = obs_dt_val.Lon.max()-obs_dt_val.Lon.min()
bdy = 0.2
box = [obs_dt_val.Lon.min()-lon_r*bdy, obs_dt_val.Lon.max()+lon_r*bdy, obs_dt_val.Lat.min()-lat_r*bdy, obs_dt_val.Lat.max()+lat_r*bdy]

dem = dem.where(dem>0)
dem = dem.where((dem.x>box[0]) & (dem.x<box[1]) & (dem.y>box[2]) & (dem.y<box[3]), drop=True)

cmap_elev = mpl.colors.ListedColormap(sns.color_palette("Greys", 10))
cmap_network = {'bcRiverForecastCenter': 'r', 'environmentAlberta': 'k', 'CRHO': 'b', 'ABE_AGG_HIST': 'g'}
cmap_dict = {'snowdepthavg':mpl.colors.ListedColormap(sns.color_palette("Blues", 12)),
             'rmse': mpl.colors.ListedColormap(sns.color_palette("Blues", 12)),
             'bias': mpl.colors.ListedColormap(sns.color_palette("RdBu_r", 12))}

# sns.set_style('whitegrid')
sns.set_context("talk", font_scale=2, rc={"lines.linewidth": 2})

#############################################################
# Make spatial plots
#############################################################

'''Above we have subset OBS and MOD to have the same stations, variables, time start/stop, and missing periods.
Below, we can caculate individual metrics, and call a common plotting script to make maps of point errors'''

# Example: Annual Average Bias
for cvar in ['t','rh','U_2m_above_srf','p','ilwr','iswr']:
    ctype = 'bias'
    da_metric = chmF.calc_bias(obs_dt_val, mod_dt_val, cvar)
    # Make plot
    cf = chmF.plot_point_metric(dem, da_metric, plot_key[cvar], ylabel_unit[cvar], cmap_dict[ctype], ctype)
    # Save Figure
    file_out = os.path.join(fig_dir, cvar+'_'+ctype+'.png')
    chmF.save_figure(cf, file_out, fig_res)

# Example: Hourly RMSE
for cvar in ['t', 'rh', 'U_2m_above_srf', 'p', 'ilwr', 'iswr']:
    ctype = 'rmse'
    da_metric = chmF.calc_rmse(obs_dt_val, mod_dt_val, cvar)
    # Make plot
    cf = chmF.plot_point_metric(dem, da_metric, plot_key[cvar], ylabel_unit[cvar], cmap_dict[ctype], ctype)
    # Save Figure
    file_out = os.path.join(fig_dir, cvar+'_'+ctype+'.png')
    chmF.save_figure(cf, file_out, fig_res)

# plt.show()



# In[11]:
#
# fig = plt.figure(figsize=(20, 20))
# ax1 = plt.axes(projection=ccrs.AlbersEqualArea())
# # ax1.set_extent(box)
# ax1.imshow(np.flipud(dem.values), extent=[np.min(dem.x), np.max(dem.x),
#                                           np.min(dem.y), np.max(dem.y)], aspect=ax1.get_aspect())
# # ax1.set_title('Elevation')
# for c_net in set(obs_dt_val.network.values):
#     lat_pts = obs_dt_val.Lat.sel(station=(obs_dt_val.where(obs_dt_val.network == c_net, drop=True).network).station).values
#     lon_pts = obs_dt_val.Lon.sel(station=(obs_dt_val.where(obs_dt_val.network == c_net, drop=True).network).station).values
#     I_not_nan = ~np.isnan(lat_pts) & ~np.isnan(lon_pts)
#     lat_pts = lat_pts[I_not_nan]
#     lon_pts = lon_pts[I_not_nan]
#
#     ax1.scatter(lon_pts, lat_pts, transform=ccrs.AlbersEqualArea(), s=50, c=cmap_network[c_net], zorder=100,
#                 label=c_net)  # yc, xc -- lists or numpy arrays
#
# # Snow Courses
# lat_pts = EC_data.Lat.values
# lon_pts = EC_data.Lon.values
# I_not_nan = ~np.isnan(lat_pts) & ~np.isnan(lon_pts)
# lat_pts = lat_pts[I_not_nan]
# lon_pts = lon_pts[I_not_nan]
#
# ax1.scatter(lon_pts, lat_pts, transform=ccrs.AlbersEqualArea(), marker='o', s=50, c='m', zorder=200,
#             label='EC Snow Course')  # yc, xc -- lists or numpy arrays
#
# ax1.add_geometries(p_sh, ccrs.AlbersEqualArea(),
#                    edgecolor='black', facecolor='none', alpha=0.5)
# plt.legend()
#
# plt.show()