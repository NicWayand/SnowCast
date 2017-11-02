import numpy as np
import xarray as xr
import sys
import os
import imp
import seaborn as sns
sns.set_context("talk",font_scale=1.5)
sns.set_style('whitegrid')

# Load in config file
#######  load user configurable paramters here    #######
# Check user defined configuraiton file
if len(sys.argv) == 1:
    sys.error('Requires one argument [configuration file]')

# Get name of configuration file/module
configfile = sys.argv[-1]

# Load in configuration file as module
X = imp.load_source('',configfile)

# Assign to local variables
data_dir = X.data_dir

AB_recent    = os.path.join(data_dir,'AB_recent','netcdf','AB_NRT.nc')
AB_POR       = os.path.join(data_dir,'AB_POR','netcdf','AB_SWE_SD_NRT.nc')
BC_CWY       = os.path.join(data_dir,'BC_NRT','netcdf','BC_NRT.nc') # Current Water Year (CWR)
BC_HIST      = os.path.join(data_dir,'BC_HIST','netcdf','BC_HIST.nc')
ABE_AGG_HIST = os.path.join(data_dir,'ABE_AGG_HIST','netcdf','ABE_AGG_HIST.nc') # MST
CRHO         = os.path.join(data_dir,'CRHO_HIST','netcdf','CRHO_1hour.nc') # MST

# In[ ]:

CRHO_time_zone = 7 # MST to UTC
ABE_AGG_HIST_time_zone = 7 # MST to UTC


# In[ ]:

merged_dir = os.path.join(data_dir,'merged')
# Make if does not exist
if not os.path.exists(merged_dir):
    os.makedirs(merged_dir)
netcdf_file_out = os.path.join(merged_dir, 'Hourly_Merged.nc')


# # Merge Networks

# In[ ]:

ds_AB_recent = xr.open_dataset(AB_recent)
ds_AB_POR    = xr.open_dataset(AB_POR)
ds_BC_CWR    = xr.open_dataset(BC_CWY)
ds_BC_HIST   = xr.open_dataset(BC_HIST)
ds_CRHO      = xr.open_dataset(CRHO)
ds_ABE_AGG   = xr.open_dataset(ABE_AGG_HIST)


# In[ ]:




# In[ ]:

# for csta in ds_ABE_AGG.staID:
#     x = ds_ABE_AGG.WindSpeed.sel(staID=csta)
#     if x.notnull().sum(dim='Time_MST')>0:
#         plt.plot(ds_ABE_AGG.Time_MST, x)
#         plt.ylim([0,360])


# In[ ]:

# Make all datasets same time zone and variable names
# All to MST
ds_CRHO['time_hrly'] = ds_CRHO['time_hrly'] + np.timedelta64(CRHO_time_zone,'h')
ds_CRHO.rename({'station':'staID','time_hrly':'Time_UTC'}, inplace=True);

ds_ABE_AGG['Time_MST'] = ds_ABE_AGG['Time_MST'] + np.timedelta64(ABE_AGG_HIST_time_zone,'h')
ds_ABE_AGG.rename({'Time_MST':'Time_UTC'}, inplace=True);


# In[ ]:

list(ds_CRHO.data_vars)


# In[ ]:

# Drop CRHO variables we are not interested in (still in orig file if we want them)
ds_CRHO = ds_CRHO[['WindDirectionatA',
 'TotalPressureUnadjustedA',
 'DownwardSolarRadiation',
 'ScalarWindSpeedA',
 'SnowDepthA',
 'AirtemperatureA',
 'IncrementalPrecipitationA',
 'UpwardSolarRadiation',
 'AirMoistureContentA',
 'UpwardTerrestrialRad',
 'DownwardTerrestrialRad',
 'SnowWaterEquivelentA']]


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:

orig_coords = ['Lat','Lon','Elevation','network','station_name'] # combine_first() can't handel coords merging


# In[ ]:

ds_AB_recent.reset_coords(orig_coords, inplace=True);
ds_AB_POR.reset_coords(orig_coords, inplace=True);
ds_BC_CWR.reset_coords(orig_coords, inplace=True);
ds_BC_HIST.reset_coords(orig_coords, inplace=True);
ds_CRHO.reset_coords(orig_coords, inplace=True);
ds_ABE_AGG.reset_coords(orig_coords, inplace=True);


# In[ ]:

# Merge AB data together (hist to midnight, and last few days)
###### This step drops coords info so leaving out most recent data for now
ds_AB_mrg = ds_AB_POR.combine_first(ds_AB_recent)
# ds_AB_mrg = ds_AB_POR


# In[ ]:

# Merge BC stations
ds_BC_merg = xr.merge([ds_BC_HIST, ds_BC_CWR])
ds_BC_HIST = None
ds_BC_CWR = None


# In[ ]:

# Merge AB and BC pillows
ds_BC_AB = xr.merge([ds_AB_mrg,ds_BC_merg])
ds_AB_mrg = None
ds_BC_merg = None


# In[ ]:

ds_BC_AB


# In[ ]:

# Rename BC AB naming to match CRHO name format
var_dict = {'AirTemperature':'AirtemperatureA','Precipitation':'CummulativePrecipitationA',
            'SWE':'SnowWaterEquivelentA','Snowdepth':'SnowDepthA'}

ds_BC_AB.rename(var_dict, inplace=True);


# In[ ]:




# In[ ]:

# Adjust units to be metric standard
ds_BC_AB['SnowWaterEquivelentA'] = ds_BC_AB.SnowWaterEquivelentA / 1000 # mm to m
ds_BC_AB['CummulativePrecipitationA'] = ds_BC_AB.CummulativePrecipitationA / 1000 # mm to m
ds_BC_AB['SnowDepthA'] = ds_BC_AB.SnowDepthA / 100 # cm to m


# In[ ]:




# In[ ]:

# Merge CRHO with other AB and BC data
ds_merged = xr.merge([ds_CRHO,ds_BC_AB])
ds_CRHO = None
ds_BC_AB = None


# In[ ]:




# In[ ]:

# for cs in ds_ABE_AGG.staID:
#     x = ds_ABE_AGG.Precipitation.sel(staID=cs)
#     if x.notnull().sum()>0:
#         plt.plot(x.Time_UTC, x.cumsum(dim='Time_UTC'), label=(str(cs)));


# In[ ]:




# In[ ]:

# Rename ds_ABE_AGG naming to match CRHO name format
var_dict = {'RealtiveHumidity':'AirMoistureContentA','WindSpeed':'ScalarWindSpeedA','WindDirection':'WindDirectionatA',
            'AirTemperature':'AirtemperatureA','Precipitation':'IncrementalPrecipitationA',}
ds_ABE_AGG.rename(var_dict, inplace=True);


# In[ ]:

# Adjust units to be metric standard
ds_ABE_AGG['IncrementalPrecipitationA'] = ds_ABE_AGG.IncrementalPrecipitationA / 1000 # mm to m


# In[ ]:

def find_common_staID(ds1,ds2):
    same = list(set(list(ds1.staID.values)).intersection(list(ds2.staID.values)))
    return(ds1.sel(staID=same).station_name)
    # Usage find_common_staID(ds_ABE_AGG,ds_merged)


# In[ ]:

ds_ABE_AGG


# In[ ]:

ds_merged


# In[ ]:




# In[ ]:




# In[ ]:

# Merge ds_merged with ds_ABE_AGG 
# Use combine_first because there are duplicate stations (met vars and snow vars)
# ds_merged_2 = xr.merge([ds_merged, ds_ABE_AGG])
ds_merged_2 = ds_ABE_AGG.combine_first(ds_merged)
ds_merged = None
ds_ABE_AGG = None


# In[ ]:




# In[ ]:

ds_merged_2


# In[ ]:




# In[ ]:




# In[ ]:

# Set coords
ds_merged_2.set_coords(orig_coords, inplace=True);


# In[ ]:

# ds_merged.Time_UTC.diff(dim='Time_UTC').plot()


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:

# Save as netcdf file
ds_merged_2.to_netcdf(netcdf_file_out)
print(netcdf_file_out)


# In[ ]:

# set(ds_merged.network.values)


# In[ ]:



