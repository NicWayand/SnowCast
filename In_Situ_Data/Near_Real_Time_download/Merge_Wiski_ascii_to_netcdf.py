import pandas as pd
import numpy as np
import xarray as xr
import sys
import os
import imp
import wget
import seaborn as sns
import glob


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

# Data network
network = 'CRHO_NRT'
download_dir = os.path.join(data_dir,network,'current')

# Netcdf file to save to
netcdf_dir   = os.path.join(data_dir,network,'netcdf')
# Make if does not exist
if not os.path.exists(netcdf_dir):
    os.makedirs(netcdf_dir)
netcdf_file_out =  os.path.join(netcdf_dir,'CRHO_NRT.nc')

# Ascii files are downloaded as Station_name__Variable_name
os.chdir(download_dir)
all_files = glob.glob('*.csv')

# Get station and variable names
all_sta = []
all_var = []
for cf in all_files:
    fname = cf.split('.')[0]
    sta_name, var_name = fname.split('__')
    all_sta.append(sta_name)
    all_var.append(var_name)
# Get unique
all_sta = list(set(all_sta))
all_var = list(set(all_var))

# Import files
ds_list = []
for csta in all_sta:

    da_list = [] # list to store all variables for this station

    for cvar in all_var:

        cf=csta+'__'+cvar+'.csv'
        if not os.path.isfile(os.path.join(download_dir,cf)):
            continue
        print(cf)

        # Load in to python
        # 2017-09-30 17:00:00
        dateparse = lambda x: pd.datetime.strptime(x, '%Y-%m-%d %H:%M:%S')
        try: # Because sometimes the file exists, but is empty
            df = pd.read_csv(cf, index_col=1, engine='python',
                             parse_dates=True, date_parser=dateparse, na_values=['NA','no value'])
            if df.shape[0]==0:
                continue
        except:
            print("Could not parse "+cf+". Something wrong with file format.")
            continue
        df.index.names = ['Time_MST']
        s = df.ix[:, 1] # Select variable (becomes a Series)
        s.name = cvar
        da = xr.DataArray.from_series(s)
        # Store
        da_list.append(da)

    # Merge into dataset
    ds = xr.merge(da_list)
    ds.coords['station'] = csta

    # Fill in Dummy variables
    for tvar in all_var:
        if not tvar in ds.data_vars:
            # Fill it
            ds[tvar] = xr.DataArray(ds[ds.data_vars.keys()[0]]*np.NaN,
                                    coords={'Time_MST':ds.Time_MST}, dims=('Time_MST'))

    # Store in list
    ds_list.append(ds)

# Concat by stations
ds_all = xr.concat(ds_list,dim='station')

# To netcdf
ds_all.to_netcdf(netcdf_file_out)
#
# import matplotlib.pyplot as plt
# for cvar in ds_all.data_vars:
#     plt.figure()
#     for csta in ds_all.station:
#         if ds_all[cvar].sel(station=csta).notnull().sum()>0:
#             plt.plot(ds_all.Time_MST, ds_all[cvar].sel(station=csta))



    # # Rename columns
    # df = df.rename(columns = AB_2_BC_var_dict)
    #
    # # df to ds
    # ds_c = xr.Dataset.from_dataframe(df)
    # # Add as coord
    # ds_c['staID'] = csta
    # ds_c.set_coords('staID',inplace=True)
    #
    # # Store as dict (if we have any data)
    # if ds_c.Time_MST.size>0:
    #     ds_list.append(ds_c)
#                 sta_list_used.append(csta)
#
# # Concat all stations for one variable into netcdf
# ds_var = xr.concat(ds_list,'staID')
# ds_var_list.append(ds_var)

