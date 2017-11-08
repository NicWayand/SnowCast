import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import sys
import os
import imp
import seaborn as sns
sns.set_context("talk",font_scale=1.5)
sns.set_style('whitegrid')

#TODO: Add constant check

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

# # Create paths

# In[ ]:

hourly_merged = os.path.join(data_dir,'merged','Hourly_Merged.nc')


# In[ ]:

QC_dir = os.path.join(data_dir,'QC')
# Make if does not exist
if not os.path.exists(QC_dir):
    os.makedirs(QC_dir)
netcdf_file_out = os.path.join(QC_dir, 'Hourly_QC.nc')


# # QC merged data

# In[ ]:

ds = xr.open_dataset(hourly_merged) #, chunks={'Time_UTC':1, 'staID':10})


# In[ ]:

#ds


# In[ ]:

#ds.WindDirectionatA.median()


# In[ ]:

# # for v in ds.data_vars:
# #     plt.figure()
# #     plt.plot(ds.Time_UTC, ds[v].values)
#plt.figure()
#plt.plot(ds.Time_UTC,ds.SnowWaterEquivelentA.values);
#plt.figure()
#plt.plot(ds.Time_UTC,ds.SnowDepthA.values);
#plt.figure()
#plt.plot(ds.Time_UTC,ds.AirtemperatureA.values);
#plt.figure()
#plt.plot(ds.Time_UTC,ds.CummulativePrecipitationA.values);
#plt.figure()
#plt.plot(ds.Time_UTC,ds.IncrementalPrecipitationA.values);


# In[ ]:




# # Quality control data here

# In[ ]:

def QC_min_max(da, vmin, vmax):
    return da.where((da <= vmax) &  (da >= vmin))


# In[ ]:

def QC_ROC(da, ROC_thress, ROC_window):
    return da.groupby('staID').apply(lambda x: remove_outliers_via_filter(x,ROC_thress,ROC_window))


# In[ ]:

def remove_outliers_via_filter(x,threshold,window):
    if ((sum(np.isnan(x.values))/len(x.values)>0.9) | np.isnan(threshold)): # Mostly nan, just return x, otherwise filter fails
        # Threshold of np.NaN indicates this method is not suitable for this type of variable (i.e. incremntal precipitation)
        #print('Less than 10% data for ',str(x.staID),str(x.name))
        return x
    else: # Have some data, apply the median filter and remove differences greater than threshold
        # Apply median filter
        temp = x.to_series().rolling(window=window, center=True).median().fillna(method='bfill').fillna(method='ffill')
        # Take difference between filter and orig data
        difference = np.abs(x.to_series() - temp)
        # Find those data values that the diff was less than or equal to the user supplied threshold
        return x.where(difference <= threshold)


# In[ ]:

## Quality Control (Hourly)

## Max and min 
min_limits = {'WindDirectionatA':0, 'ScalarWindSpeedA':0, 'AirMoistureContentA':5,'SnowWaterEquivelentA':0, 'SnowDepthA':0, 'CummulativePrecipitationA':0, 'IncrementalPrecipitationA':0, 'AirtemperatureA':-40}
max_limits = {'WindDirectionatA':360, 'ScalarWindSpeedA':30, 'AirMoistureContentA':100,'SnowWaterEquivelentA':3, 'SnowDepthA':5.5, 'CummulativePrecipitationA':3, 'IncrementalPrecipitationA':60/1000, 'AirtemperatureA':50}

for cvar in min_limits.keys():
    print('QC-min-max: '+str(cvar))
    ds[cvar] = QC_min_max(ds[cvar], min_limits[cvar], max_limits[cvar]) 


# In[ ]:

## ROC - Use median filter to find values
ROC_thress = {'AirMoistureContentA':60,
              'SnowWaterEquivelentA':0.2, 
              'SnowDepthA':0.5, 
              #'CummulativePrecipitationA':20/1000, 
              'AirtemperatureA':10} # (unit/hr)
ROC_window = {'AirMoistureContentA':6,
              'SnowWaterEquivelentA':4, 
              'SnowDepthA':4, 
              #'CummulativePrecipitationA':48, 
              'AirtemperatureA':6} # dt (hrs)

for cvar in ROC_thress.keys():
    print('QC-ROC: '+str(cvar))
    ds[cvar] = QC_ROC(ds[cvar], ROC_thress[cvar], ROC_window[cvar]) 


# In[ ]:

# import mpld3
# mpld3.enable_notebook()


# In[ ]:

# # Aggregate to daily using mean of cummulated (to help remove noise)

# # Precip
# Daily_cum_precip = ds['CummulativePrecipitationA'].resample(freq='D',dim='Time_UTC',how='median',label='left')
# ds['Daily_QC_Cumulative_Precipitation'] = Daily_cum_precip.rename({'Time_UTC':'Daily_Time_UTC'})

# # SWE
# Daily_SWE = ds['QC_Snow_Water_Equivalent'].resample(freq='D',dim='Time_UTC',how='mean',label='left')
# ds['Daily_QC_Snow_Water_Equivalent'] = Daily_SWE.rename({'Time_UTC':'Daily_Time_UTC'})


# In[ ]:




# In[ ]:

#DS = '2012-10'
#DE = '2013-09'


# In[ ]:

#plt.plot(ds.Time_UTC.sel(Time_UTC=slice(DS,DE)),ds.SnowWaterEquivelentA.sel(Time_UTC=slice(DS,DE)).values);


# In[ ]:

#ds.AirMoistureContentA.plot.hist(bins=1000);


# In[ ]:

#plt.plot(ds.Time_UTC.sel(Time_UTC=slice(DS,DE)),ds.AirMoistureContentA.sel(Time_UTC=slice(DS,DE)).values);


# In[ ]:

#ds.ScalarWindSpeedA.plot.hist(bins=1000);


# In[ ]:

#plt.plot(ds.Time_UTC.sel(Time_UTC=slice(DS,DE)),ds.ScalarWindSpeedA.sel(Time_UTC=slice(DS,DE)).values);


# In[ ]:

#ds.WindDirectionatA.plot.hist(bins=1000);


# In[ ]:

#plt.plot(ds.Time_UTC.sel(Time_UTC=slice(DS,DE)),ds.WindDirectionatA.sel(Time_UTC=slice(DS,DE)).values);


# In[ ]:

#ds.SnowDepthA.plot.hist(bins=500);


# In[ ]:

#plt.plot(ds.Time_UTC.sel(Time_UTC=slice(DS,DE)),ds.SnowDepthA.sel(Time_UTC=slice(DS,DE)).values);


# In[ ]:

#plt.plot(ds.Time_UTC.sel(Time_UTC=slice(DS,DE)),ds.AirtemperatureA.sel(Time_UTC=slice(DS,DE)).values);


# In[ ]:




# In[ ]:

#plt.plot(ds.Time_UTC.sel(Time_UTC=slice(DS,DE)),ds.CummulativePrecipitationA.sel(Time_UTC=slice(DS,DE)).values);


# In[ ]:

#(ds.IncrementalPrecipitationA.where(ds.IncrementalPrecipitationA>0,drop=True)*1000).plot.hist(bins=100);


# In[ ]:

#plt.plot(ds.Time_UTC, ds.IncrementalPrecipitationA.cumsum(dim='Time_UTC').values*1000);
# plt.plot(ds.Time_UTC.sel(Time_UTC=slice(DS,DE)), ds.IncrementalPrecipitationA.sel(Time_UTC=slice(DS,DE)).values*1000);


# In[ ]:

#plt.plot(ds.Time_UTC.sel(Time_UTC=slice(DS,DE)), ds.IncrementalPrecipitationA.sel(Time_UTC=slice(DS,DE)).cumsum(dim='Time_UTC').values);


# In[ ]:

# Save as netcdf file
ds.to_netcdf(netcdf_file_out)
print(netcdf_file_out)


# In[ ]:



