#!/bin/bash
OMP_NUM_THREADS=24

# Run info
cd /home/nwayand/snow_models/output_CHM/SnowCast/CHM_Configs

#./CHM -f GEM_west_fortress.json #-c config.option.startdate:$start_date 
./CHM -f HRDPS_CHM_Historical.json

# Run plotting scripts
echo "Running Plotting scripts"
/home/nwayand/custom/anaconda2/bin/python /home/nwayand/SnowCast/Post_Processing/plot_scripts/Plot_SnowCast.py ../Path_Config.py HRDPS_CHM_Historical

echo DONE!!!


