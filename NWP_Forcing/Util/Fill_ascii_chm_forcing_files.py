import xarray as xr
import os
import glob
import imp
import sys
import numpy as np
import pandas as pd
import datetime
import json
import time
import utm
import matplotlib.pyplot as plt
###
# Experimental cache option to speed up dask calls
# import cachey
# from dask.cache import Cache
# cache = Cache(10e9)
# cache.register()
###
start_time = time.time()
# Hack to force datetimes to display in GMT/UTC (numpy 1.11.1 has fixed this but other dependent modules (pynio) can't handel numpy 1.11.1)
os.environ['TZ'] = 'GMT'
time.tzset()

# Load in config file
#######  load user configurable paramters here    #######
# Check user defined configuraiton file
if len(sys.argv) == 1:
    sys.error('requires one argument [configuration file] (i.e. python GRIB2_to_CHM_forcing.py forcing_config.py')

# Get name of configuration file/module
configfile = sys.argv[-1]

# Load in configuration file as module
X = imp.load_source('',configfile)

output_dt = 1 # Hours

# Assign to local variables
netcdf_dir = X.netcdf_dir
ascii_dir   = X.ascii_dir
Forcing_config_file = X.Forcing_config_file
lat_r = X.lat_r
lon_r = X.lon_r
var_dic = X.var_dic
local_time_offset = X.local_time_offset
coordsystem = X.coordsystem
if coordsystem=='pro':
    utm_zone = X.utm_zone

# Output dir
output_dir = ascii_dir + '_filled'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Move to input
os.chdir(ascii_dir)

all_files = sorted(set(glob.glob('*.chm')))

class ascii_file(object):


    def __init__(self, cfile=None):
        # Load it in
        cfile = open(cfile,'r')
        df = pd.read_csv(cfile, sep="\t", parse_dates=True)
        cfile.close()
        df.set_index('datetime', inplace=True)
        df.index = pd.to_datetime(df.index)
        self.df = df

    def add_missing_timesteps(self, freq='H'):
        # Reindex to have continuous time steps
        self.df_c = self.df.reindex(pd.date_range(self.df.index[0], self.df.index[-1],
                                          freq=freq))

    def fill_missing_variables(self, method='linear'):
        # Simple linear interpolation along time
        if method == 'linear' or method == 'spline':
            self.df_c = self.df_c.interpolate(method=method, axis=0).ffill().bfill()
        else:
            raise ValueError('Method not found.')

        assert not self.df_c.isnull().values.any()

    def write_to_ascii(self, cfile=None):
        self.df_c.index.name = 'datetime' # This gets dropped so add back in
        file_out = open(os.path.join(output_dir, cfile), 'w')
        self.df_c.to_csv(file_out, sep='\t', date_format='%Y%m%dT%H%M%S')
        file_out.close()

for cf in all_files:
    print(cf)
    cO = ascii_file(cf)
    cO.add_missing_timesteps(freq='H')
    cO.fill_missing_variables(method='linear')
    cO.write_to_ascii(cf)

    # # Test Plot
    # plt.figure()
    # plt.plot(cO.df.index, cO.df['t'],'k*',linewidth=6)
    # plt.plot(cO.df_c.index, cO.df_c['t'], 'r')
    # plt.plot(cO.df_c.index, cO.df_c['t'], 'ob')
    #
    # plt.show()



