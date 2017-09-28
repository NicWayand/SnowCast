#!/bin/bash

rsync -asv /home/nwayand/SnowCast/CHM_Configs/forecast_CRHO_spinup/figures/* root@www.snowcast.ca:/var/www/html/static/HRDPS/
rsync -asv /home/nwayand/SnowCast/CHM_Configs/GDPS_Current/figures/* root@www.snowcast.ca:/var/www/html/static/GDPS/

#sshpass -p PASSWORDHERE rsync /home/nwayand/SnowCast/CHM_Configs/forecast_CRHO_spinup/figures/* root@www.snowcast.ca:/var/www/html/static/HRDPS/


