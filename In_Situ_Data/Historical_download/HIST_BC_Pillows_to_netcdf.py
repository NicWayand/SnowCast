import numpy as np
import pandas as pd
import xarray as xr
import sys
import os
import imp
import wget
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
git_dir   = X.git_dir

# # Create paths

# In[3]:

# Data network
network = 'BC_HIST'

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
netcdf_file_out =  os.path.join(netcdf_dir,'BC_HIST.nc')

# Metadata for AB pillows 
meta_file         = 'BC_Station_Metadata.csv'
meta_file_path    = os.path.join(git_dir,'In_Situ_Data','metadata',meta_file)

# # Download Historical BC data (Updated at end of water year)

# In[4]:

os.chdir(download_dir)
# File format
# 
# File names change only by last water year, so figure this out based on current date
l_WY = '2016'
BC_files = ['Raw Snow Water Equivalent (mm) October 1 2011 to September 30 '+l_WY+'.csv', 
           'Raw Cumulative Precipitation (mm) October 1 2011 to September 30 '+l_WY+'.csv',
            'Raw Snow Depth (cm) October 1 2011 to September 30 '+l_WY+'.csv',
           'Raw Air Temperature (Degrees Celsius) October 1 2011 to September 30 '+l_WY+'.csv']
Var_names = ['SWE','Precipitation','Snowdepth','AirTemperature']
Var_units = ['mm','mm','cm','C']
c_network = 'bcRiverForecastCenter'


# In[5]:

# Remove previous files
for cfile in BC_files:
    try:
        os.remove(cfile)
    except OSError:
        pass


# In[6]:

# Download newest files
url_base         = 'http://bcrfc.env.gov.bc.ca/data/asp/'
[wget.download(url_base+cfile) for cfile in BC_files]


# In[7]:

# Import metadata for each station
metadata = pd.read_csv(meta_file_path,index_col=1,delimiter=',',encoding='utf-8')


# In[8]:

# Import each data variable
ds_dict   = {}
unit_dict = {}
for (i,cf) in enumerate(BC_files):
    print(cf)
    
    # Load in to python
    df = pd.read_csv(cf,index_col=0, skipfooter=1, engine='python', parse_dates=True)
    var_name_full  = Var_names[i]
    var_units = Var_units[i]
    df.index.names = ['Time_UTC']
    
#     # Check for error in PC.csv
#     if cf=='PC.csv':
#         print('fixing Muskwa-Kechik error')
#         df = df.rename(columns = {'4A34P Muskwa-Kechika':'4A34P Dowling Creek'})

    # Store as dict
    ds_dict[var_name_full] = df
    unit_dict[var_name_full] = var_units


# In[9]:

# Merge into netcdf
ds = xr.Dataset(ds_dict)
ds.rename({'dim_1':'staID'},inplace=True) # rename time
ds['staID']        = [str(x).split(' ')[0] for x in ds.staID.values]


# In[10]:

## ADD UNITS
# Add variable attributes (units), and fix variable names (remove spaces)
for cvar in ds.data_vars:
    # add units as attributes
    ds.get(cvar).attrs['unit']   = unit_dict[cvar]


# In[11]:

## Add station metadata
ds['station_name'] = xr.DataArray(metadata['station'],coords={'staID':metadata.index}, dims='staID')
ds['Lat'] = xr.DataArray(metadata['latitude'],coords={'staID':metadata.index}, dims='staID')
ds['Lon'] = xr.DataArray(metadata['longitude'],coords={'staID':metadata.index}, dims='staID')
ds['Elevation'] = xr.DataArray(metadata['elevation'],coords={'staID':metadata.index}, dims='staID')


# In[12]:

# Move to coords
ds.set_coords(['station_name','Lat','Lon','Elevation'], inplace=True)
ds


# In[ ]:




# In[13]:

# Reindex to have continous time steps
Time_UTC_new = np.arange(ds.Time_UTC.values[0], ds.Time_UTC.values[-1], dtype='datetime64[h]')
ds_fill = ds.reindex({'Time_UTC':Time_UTC_new})


# In[ ]:




# In[ ]:




# In[14]:

# Add Network
ds_fill.coords['network'] = xr.DataArray([c_network for x in ds_fill.staID], dims='staID')


# In[15]:

# Save as netcdf file
ds_fill.to_netcdf(netcdf_file_out)


# In[ ]:



