#!/bin/bash
#MP_NUM_THREADS=4

# Caller script to run CHM for all GEM hours

# Run info
cd /home/nwayand/snow_models/output_CHM/SnowCast

#./CHM -f GEM_west_fortress.json #-c config.option.startdate:$start_date 
#./CHM -f GEM_west_CRHO.json

# Run plotting scripts
echo "Running Plotting scripts"
/home/nwayand/custom/anaconda2/bin/python ./plot_scripts/Plot_SnowCast.py
echo "Plotting complete"

# Upload to server
# TODO: build into python plotting scripts like this (http://stackoverflow.com/questions/12613797/python-script-uploading-files-via-ftp)
echo "Uploading to server"
#cd /home/nwayand/snow_models/output_CHM/SnowCast/forecast_spinup/figures
cd /home/nwayand/snow_models/output_CHM/SnowCast/forecast_CRHO_spinup/figures
../../ftp_upload.sh density.png files
../../ftp_upload.sh snowdepth_48h.png files
../../ftp_upload.sh 48_snowdepth_change.png files
echo "Upload complete"

# Convert vtus to geotifs
/home/nwayand/custom/anaconda2/bin/python /home/nwayand/snow_models/CHM/tools/vtu2geo/main.py vtu2geo_config.py

# Convert vtus to KMZ and uploa
# Convert vtus to KMZ and upload
cd /home/nwayand/snow_models/output_CHM/SnowCast/KML_files
./batch_tif_2_KML_UPLOAD.sh /home/nwayand/snow_models/output_CHM/SnowCast/KML_files snowdepthavg

echo DONE!!!


