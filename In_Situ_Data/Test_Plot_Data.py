import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import sys
import imp
import os
import seaborn as sns
plt.rcParams.update({'figure.max_open_warning': 0})


ds_m = xr.open_dataset(r'/media/data2/SnowCast_station_data/merged/Hourly_Merged.nc')
ds_q = xr.open_dataset(r'/media/data2/SnowCast_station_data/QC/Hourly_QC.nc')

# General plotting settings
sns.set_style('ticks')
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})

# time
#ds_m = ds_m.sel(Time_UTC=slice("2014-10-01","2017-09-28"))
#ds_q = ds_q.sel(Time_UTC=slice("2014-10-01","2017-09-28"))

# vars_to_plot = ['SnowWaterEquivelentA','SnowDepthA','AirtemperatureA']
vars_to_plot = ['IncrementalPrecipitationA']
# network = 'ABE_AGG'
network = 'CRHO'

ds_m = ds_m[vars_to_plot]
ds_q = ds_q[vars_to_plot]

if network!='all':
    ds_m = ds_m.where(ds_m.network==network, drop=True)
    ds_q = ds_q.where(ds_q.network==network, drop=True)

for cvar in vars_to_plot:

    plt.figure()
    plt.title('Merged')
    #plt.plot(ds_m.Time_UTC, ds_m[cvar].sel(staID=['05BJ805' ,'05CA805' ,'2A32P' ,'2C14P' ,'PWL']))
    for csta in ds_m.staID:
        plt.plot(ds_m.Time_UTC, ds_m[cvar].sel(staID=csta),
                 label=ds_m.sel(staID=csta).station_name.values)
    plt.legend()

    plt.figure()
    plt.title('QC')
    #plt.plot(ds_q.Time_UTC, ds_q[cvar].sel(staID=['05BJ805' ,'05CA805' ,'2A32P' ,'2C14P' ,'PWL']))
    for csta in ds_q.staID:
        plt.plot(ds_q.Time_UTC, ds_q[cvar].sel(staID=csta),
                 label=ds_q.sel(staID=csta).station_name.values)
    plt.legend()

plt.show()


