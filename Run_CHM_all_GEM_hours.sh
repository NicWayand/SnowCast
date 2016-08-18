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


