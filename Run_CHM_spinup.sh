#!/bin/bash
#MP_NUM_THREADS=4

# Caller script to run CHM for all GEM hours

# Run info
cd /home/nwayand/snow_models/output_CHM/SnowCast

#./CHM -f GEM_west_fortress.json #-c config.option.startdate:$start_date 
./CHM -f GEM_west_CRHO.json

# Run plotting scripts
echo "Running Plotting scripts"
/home/nwayand/custom/anaconda2/bin/python ./plot_scripts/Plot_SnowCast.py
echo "Plotting complete"

# Upload to server
# TODO: build into python plotting scripts like this (http://stackoverflow.com/questions/12613797/python-script-uploading-files-via-ftp)
echo "Uploading to server"
#cd /home/nwayand/snow_models/output_CHM/SnowCast/forecast_spinup/figures
cd /home/nwayand/snow_models/output_CHM/SnowCast/forecast_CRHO_spinup/figures
./ftp_upload.sh 48_snowfall.png
./ftp_upload.sh snowdepth_48h.png
./ftp_upload.sh 48_snowdepth_change.png
echo "Upload complete"


