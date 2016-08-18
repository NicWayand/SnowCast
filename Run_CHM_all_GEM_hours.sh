#!/bin/bash
MP_NUM_THREADS=4

# Caller script to run CHM for all GEM hours

# Run info
cd /home/nwayand/snow_models/output_CHM/SnowCast
#start_date = "20160817T200000"
#end_date   = ""

# Move to main dir
#cd $main_dir
#echo pwd

./CHM -f GEM.json #-c config.option.startdate:$start_date 

# Convert output vtu files to tif files
python /home/nwayand/snow_models/CHM/tools/vtu2geo/main.py vtu2geo_config.py

