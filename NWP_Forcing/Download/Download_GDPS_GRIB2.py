import os
import numpy as np
import wget
import sys
import imp
import time

# Load in config file
#######  load user configurable paramters here    #######
# Check user defined configuraiton file 
if len(sys.argv) == 1:
    sys.error('Download_RDPS_GRIB2.py requires one argument [configuration file] (i.e. python Download_HRDPS_GRIB2.py download_config.py')

# Get name of configuration file/module
configfile = sys.argv[-1]

# Load in configuration file as module
X = imp.load_source('',configfile)

# Assinge to local variables
download_dir = X.download_dir
Init_H = X.Init_H
Variable = X.Variable
domain = X.domain 

# Move to input if exists (make if not)
if not os.path.isdir(download_dir):
    os.mkdir(download_dir)
os.chdir(download_dir)

# REMOVE ALL previous files
print('Removing ALL previous grib2 files')
files = os.listdir(download_dir)
for file in files:
    os.remove(file)

# Define format of date in file names
YYYYMMDDHH = time.strftime("%Y%m%d")+Init_H

# Fixed parameters for path and file name syntax
# CMC_glb_Variable_LevelType_Level_projection_YYYYMMDDHH_Phhh.grib2
# CMC_glb_TMP_ISBL_925_latlon.24x.24_2010090800_P042.grib2
base_url = 'http://dd.weather.gc.ca/model_gem_global'
filetype = 'grib2'
lat_lon = 'lat_lon'
Forc_H   = np.arange(3, 145, 3)
sep      = '/'
office     = 'CMC'
system     = 'glb'
Phhh       = 'P'
ending     = '.grib2'
us         = '_'

# http://dd.weather.gc.ca/model_gem_global/25km/grib2/lat_lon/00/003/

# For each forecast hour (Forc_H)
for c_Forc_H in Forc_H:
    # Convert to string and pad with zeros
    FH_s = np.array_str(c_Forc_H).zfill(3)

    # Get current path
    cpath = base_url + sep + domain + sep + filetype + sep + lat_lon + sep + Init_H + sep + FH_s + sep
    #print(cpath)

    # For each variable
    for cVar in np.arange(0,len(Variable)):
        
        # Make file path
        # CMC_glb_TMP_ISBL_925_latlon.24x.24_2010090800_P042.grib2
        cfile = office + us + system + us + Variable[cVar] + us + 'latlon.24x.24' + us + YYYYMMDDHH + us + Phhh + FH_s + ending
        print(cfile)
        print(cpath+cfile)
        print(" ")
        # download file
        filename = wget.download(cpath+cfile)
