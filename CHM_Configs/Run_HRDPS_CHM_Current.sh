#!/bin/bash
OMP_NUM_THREADS=5

# Caller script to run CHM for all GEM hours
N_threads=10 # number of threads to use for vtu2geo parallel

# Run info
cd /home/nwayand/snow_models/output_CHM/SnowCast/CHM_Configs

./CHM -f HRDPS_CHM_Current.json

# Run plotting scripts
echo "Running Plotting scripts"
/home/nwayand/custom/anaconda2/bin/python /home/nwayand/SnowCast/Post_Processing/plot_scripts/Plot_SnowCast.py forecast_CRHO_spinup
echo "Plotting complete"

# Upload to server
# TODO: build into python plotting scripts like this (http://stackoverflow.com/questions/12613797/python-script-uploading-files-via-ftp)
echo "Uploading to server"
#cd /home/nwayand/snow_models/output_CHM/SnowCast/forecast_spinup/figures
cd /home/nwayand/snow_models/output_CHM/SnowCast/CHM_Configs/forecast_CRHO_spinup/figures
/home/nwayand/SnowCast/Web_Upload/ftp_upload.sh density.png files
/home/nwayand/SnowCast/Web_Upload/ftp_upload.sh snowdepth_48h.png files
/home/nwayand/SnowCast/Web_Upload/ftp_upload.sh 48_snowdepth_change.png files
echo "Upload complete"

# Convert vtus to geotifs
#cd /home/nwayand/snow_models/output_CHM/SnowCast
#/home/nwayand/custom/anaconda2/bin/python /home/nwayand/snow_models/CHM/tools/vtu2geo/create_parallel_configs.py vtu2geo_config.py
#/home/nwayand/snow_models/CHM/tools/vtu2geo/run_vtu2geo_in_parallel.sh /home/nwayand/snow_models/output_CHM/SnowCast/vtu_config $N_threads

# Convert vtus to KMZ and upload
#cd /home/nwayand/snow_models/output_CHM/SnowCast/KML_files
#./batch_tif_2_KML_UPLOAD.sh /home/nwayand/snow_models/output_CHM/SnowCast/KML_files snowdepthavg

echo DONE!!!


