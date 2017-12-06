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
###
# Experimental cache option to speed up dask calls
import cachey 
from dask.cache import Cache
cache = Cache(10e9)
cache.register()
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

# Assign to local variables
netcdf_dir = X.netcdf_dir
ascii_dir   = X.ascii_dir
hgt_file = X.hgt_file
Forcing_config_file = X.Forcing_config_file
var_dic = X.var_dic
lat_r = X.lat_r
lon_r = X.lon_r
local_time_offset = X.local_time_offset
coordsystem = X.coordsystem
if coordsystem=='pro':
    utm_zone = X.utm_zone

# Hgt not included in GEM files so import here and merge later
#ds_hgt = xr.open_dataset(hgt_file,engine='netcdf4')
# Drop time (same for all times)
#ds_hgt = ds_hgt.isel(time=0)

# Loop by year and month (expand as needed)
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
years = ['2012','2013','2014','2015','2016','2017']

first_month = True # Flag to write header to csv files if first month

for cy in years:
    for cm in months:
	# Move to input
	os.chdir(netcdf_dir)

        file_exp = 'GEM_rockies_'+str(cy)+str(cm)+'*.nc'
        print 'Looking for... '+file_exp
        cfiles = glob.glob(file_exp)
        print "found "+str(len(cfiles))+" files"

        # Break out if no files found
	if len(cfiles)==0:
		continue

	# Sort file names (open_mfdataset _should_ handle this, but this is a hack to put in correct order)
	cfiles = sorted(cfiles)

	# Load all hourly WRF files
	ds = xr.open_mfdataset(cfiles,concat_dim='time',engine='netcdf4')

	# Merge with height
        #ds = xr.merge([ds,ds_hgt])

	# Rename variables
	ds.rename(var_dic,inplace=True)

	# Remove uneeded vars
	#ds = ds.drop(['])

	print 'Converting units to CHM units'
	##### Convert units to CHM requirements

	# longitude
	ds['longitude'] =  ds.longitude - 360

	# replace x,y (in distance (m)) to indices (wgrib does this for some reason!)
	ds['x']=np.arange(0,ds.x.size)
	ds['y']=np.arange(0,ds.y.size)

	# Apply a bias correction (optional)
	#ds['Qli'] = ds.Qli + 30
	#print "Warning, applying +30 ilwr bias correctoin"

        # Drop sigma dims
        ds['t'] = ds.t.isel(sigma=0)
        ds['rh'] = ds.rh.isel(sigma=0)
        ds['GZ'] = ds.GZ.isel(sigma=0)
        ds['u'] = ds.u.isel(sigma=0)
        ds['vw_dir'] = ds.vw_dir.isel(sigma=0)

	# Relative humidity (Based on Tetens' formula (1930))
        ds['rh'] = ds.rh*100 # fraction to %    

	# Pressure (after RH calcs)
	ds['press'] = ds['press'] / 100 # Pa to hPa

	# Precipitaiton (liquid and solid) accumulated (m) to incremental (mm)
        ds_p = ds.PR.dif(dim='time')
        # Setp values just below zero to zero
        ds_p.values[ds_p.values<0] = 0
        # First value is unknown (downside of saving as accum...) so we set it to -9999
        ds['p'] = xr.concat([ds.PR[0,:,:]*0-9999,ds_p],dim='time').transpose('time','y','x')

        # Geopotentila to height
        ds['HGT_surface'] = ds.GZ/9.81

	# Rename time
	ds.rename({'time':'datetime'},inplace=True)

	# Adjust to local time zone (i.e. from UTC to MST, local_time_offset should = -7)
	ds['datetime'] = pd.to_datetime(ds.datetime.values) + datetime.timedelta(hours=local_time_offset)

	# Move to ascii dir
	if not os.path.isdir(ascii_dir):
	    os.mkdir(ascii_dir)
	os.chdir(ascii_dir)

	# Extract grid cells we want to export
	print 'Extracting cells within lat/long box'
	ds_new = ds.where((ds.latitude>lat_r[0]) & (ds.latitude<lat_r[1]) & (ds.longitude>lon_r[0]) & (ds.longitude<lon_r[1]), drop=True)
	# Clear mem
	ds = None

	print 'Loading into memory'
 	# Load into memory
	ds_new.load()
	print 'done loading'

	# Get small dataset of lat and long
	grid = ds_new[['latitude', 'longitude','HGT_surface']]

	# Initalize dictionary of metadata
	meta = {}

	# Loop through each point within our lat long box
	i_range = range(ds_new.coords['x'].size)
	j_range = range(ds_new.coords['y'].size)
	print "Found " + str(ds_new.coords['x'].size) + "i and " + str(ds_new.coords['y'].size) + "j cells."
	for i in i_range:
	    for j in j_range:	
		sub_grid = grid.isel(x=i, y=j)
		# Check if in domain
		if ((sub_grid.latitude>lat_r[0]) & (sub_grid.latitude<lat_r[1]) & (sub_grid.longitude>lon_r[0]) & (sub_grid.longitude<lon_r[1])):	
		    print "processing point " + str(i) + " " + str(j)
		    # Select point
                    print "grabbing ds at point"
		    sub_ds = ds_new.isel(x=i, y=j)
		    print "converting to dataframe from Dataset"
		    # Convert to dataframe
		    df = sub_ds.to_dataframe()
                    # clean up dataset
                    sub_ds = None
		    # Drop unwanted vars/coords 
		    df.drop(['HGT_surface','latitude','longitude','x','y'], axis=1, inplace=True)

		    # Check if first month
		    if first_month:
			    # Open file for writing
			    of = open('point_'+str(i)+'_'+str(j)+'.chm','w')
		 	    
			    # Write to ascii CHM format with header
 	                    df.to_csv(of,sep='\t',date_format='%Y%m%dT%H%M%S')

			    # Close file and clean up dataframe
			    of.close()
                            df = None

			    # Build dictionary.json of point metadata
			    coords = []
			    if coordsystem=='geo':
				coords.append(float(sub_grid.longitude.values))
				coords.append(float(sub_grid.latitude.values))
			    elif coordsystem=='pro':
				coords = utm.from_latlon(sub_grid.latitude,sub_grid.longitude, utm_zone)
			    else:
				print 'Unknown coords = ' + coordsystem

			    meta['point_'+str(i)+'_'+str(j)] = {'file':ascii_dir+'/point_'+str(i)+'_'+str(j)+'.chm',
				"longitude":coords[0],
				"latitude":coords[1],
				"elevation":sub_grid.HGT_surface.data.tolist(),
				"filter": {
				"scale_wind_speed": {
				"variable": "u",
				"Z_F": 10
				}}}
		    else:
			    # Open file for writing
                            of = open('point_'+str(i)+'_'+str(j)+'.chm','a')


		      	    # Write to ascii CHM format without header, append
    	                    df.to_csv(of,sep='\t',date_format='%Y%m%dT%H%M%S', mode='a', header=False)
 
			    # Close file and clean up dataframe
			    of.close()
			    df = None

		    # Clean up memory
                    sub_grid = None                     

        # End of each x/y
        ds_new = None
        # Continue of current month
        if first_month:
		# Write out json file with metadata for these points
		print("Writing out Lam West forcing json file")
		with open(Forcing_config_file, 'w') as outfile:
	    		json.dump(meta, outfile,indent=4, sort_keys=True)

   	# End of current month
        first_month = False

print("--- %s minutes ---" % ((time.time() - start_time)/60))
print 'finished!'
