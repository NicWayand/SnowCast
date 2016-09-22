#!/bin/bash
#MP_NUM_THREADS=4

# Caller script to run CHM for all GEM hours

# Run info
cd /home/nwayand/snow_models/output_CHM/SnowCast

./CHM -f GEM_west_fortress.json #-c config.option.startdate:$start_date 

# Convert output vtu files to tif files
#python /home/nwayand/snow_models/CHM/tools/vtu2geo/main.py vtu2geo_config.py
echo Skipping vtu2geo, need to update!

# Run ipython plotting scripts
#runipy Plot_CHM_output.ipynb Spatial_plots.ipynb
runipy ./plot_scripts/Plot_scatter_CHM_OBS.ipynb Scatter_plots.ipynb
runipy ./plot_scripts/Plot_timeseries_at_CRHO_stations.ipynb Timeseries_plots.ipynb 

# Push ipython notebook plots to github
#http://askubuntu.com/questions/306176/hourly-git-push
./Git_commit_and_push_master.sh

