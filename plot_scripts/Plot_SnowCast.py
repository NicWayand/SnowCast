import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import time
import vtk
from vtk.util import numpy_support as VN
import matplotlib.tri as tri
import vtu_functions as vfunc
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import datetime

# Plot settings
import seaborn as sns
sns.set_style('whitegrid')
sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})

# Dir to vtu files
#vtu_dir   = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/SnowCast/forecast_spinup_highres/meshes')  
#fig_dir   = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/SnowCast/forecast_spinup_highres/figures')
#vtu_dir   = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/SnowCast/forecast_spinup/meshes')
#fig_dir   = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/SnowCast/forecast_spinup/figures')
vtu_dir   = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/SnowCast/forecast_CRHO_spinup/meshes')
fig_dir   = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/SnowCast/forecast_CRHO_spinup/figures')
prefix = 'SC'

#vtu_dir   = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/GEM_CRHO/test_2/meshes')
#fig_dir   = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/GEM_CRHO/test_2/figures')
#prefix    = 'GemK'

local_time_offset = -7

# Make fig dir
if not os.path.isdir(fig_dir):
    os.mkdir(fig_dir)

# Input
var_name = 'snowdepthavg'

# Move to vtu dir
os.chdir(vtu_dir)
# Get list of files
vtu_files = np.sort(glob.glob(prefix+'*.vtu'))
# file to plot
cfile = vtu_files[-1]
# Get time
time_stamp = vfunc.get_vtu_time([cfile],prefix)


# Adjust time zone
time_stamp = pd.to_datetime(time_stamp) + datetime.timedelta(hours=local_time_offset)

#print time_stamp
# Get mesh
cmesh = vfunc.get_mesh(cfile)
z = vfunc.get_face_var(cmesh,var_name)
# Get tri info
tri_info = vfunc.get_triangle(cmesh)

# Move to fig dir
os.chdir(fig_dir)

# Plot
#f1,ax1 = plt.subplots(ncols=1)
#f1.set_size_inches(10, 10)
#ax1.set_title('Snow depth. '+pd.to_datetime(time_stamp[0]).strftime('%Y-%m-%d %H'))
#p1 = plt.tripcolor(tri_info['X'], tri_info['Y'], tri_info['triang'],facecolors=z)
#b1 = plt.colorbar(p1)
## Save to file
#f1.savefig('snowdepth_48h.png',bbox_inches='tight',dpi=300)

# Plot setup
def make_map(projection=ccrs.PlateCarree()):
    fig, ax = plt.subplots(subplot_kw=dict(projection=projection))
    fig.set_size_inches(20, 20)
    #ax.coastlines(resolution='100m', zorder=1)
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    return fig, ax

# Plot on map
fig, ax = make_map()
#kw = dict(marker='.', linestyle='-', alpha=0.25, color='darkgray', zorder=0)
#lines = ax.triplot(tri_info['X'], tri_info['Y'], **kw)
ax.set_title('Snow depth. '+pd.to_datetime(time_stamp[0]).strftime('%Y-%m-%d %H:%M:%S MST'))
p1 = ax.tripcolor(tri_info['X'], tri_info['Y'], tri_info['triang'],facecolors=z*100,cmap='Blues',vmin=0, vmax=200)
b1 = fig.colorbar(p1)
b1.ax.set_ylabel('Snow depth (cm)')

sns.set_palette("BuGn_r")
fig.savefig('snowdepth_48h.png',bbox_inches='tight',dpi=300)

# Input
var_name = 'p_snow'

# Move to vtu dir
os.chdir(vtu_dir)
# Get list of files
vtu_files = np.sort(glob.glob(prefix+'*.vtu'))
# file to plot
cfile = vtu_files[-48:]

# Get time
time_stamp = vfunc.get_vtu_time(cfile,prefix)

# Adjust time zone
time_stamp = pd.to_datetime(time_stamp) + datetime.timedelta(hours=local_time_offset)

# Get data
df_p_snow = vfunc.get_multi_mesh_var_dask(cfile,var_name,time_stamp)

# Sum last 48 hours (mm to cm)
snowfall = df_p_snow.sum(axis=1)/100

# Move to fig dir
os.chdir(fig_dir)

# Plot on map
fig2, ax = make_map()
#kw = dict(marker='.', linestyle='-', alpha=0.25, color='darkgray', zorder=0)
#lines = ax.triplot(tri_info['X'], tri_info['Y'], **kw)
ax.set_title('Snow accumulation.\n '+pd.to_datetime(time_stamp[0]).strftime('%Y-%m-%d %H:%M:%S MST')+' to '+pd.to_datetime(time_stamp[-1]).strftime('%Y-%m-%d %H:%M:%S MST'))
p1 = ax.tripcolor(tri_info['X'], tri_info['Y'], tri_info['triang'],facecolors=snowfall,cmap='Blues')
b1 = fig2.colorbar(p1)
b1.ax.set_ylabel('liquid water equivalent (cm)')
sns.set_palette("BuGn_r")
fig2.savefig('48_snowfall.png',bbox_inches='tight',dpi=300)

## 48 accumued snow depth
# Input
var_name = 'snowdepthavg'

# Move to vtu dir
os.chdir(vtu_dir)
# Get list of files
vtu_files = np.sort(glob.glob(prefix+'*.vtu'))
# file to plot
cfile = [vtu_files[-48],vtu_files[-1]]

# Get time
time_stamp = vfunc.get_vtu_time(cfile,prefix)

# Adjust time zone
time_stamp = pd.to_datetime(time_stamp) + datetime.timedelta(hours=local_time_offset)

# Get data
df_snowdepth = vfunc.get_multi_mesh_var_dask(cfile,var_name,time_stamp)

# Take diff in 48 hrs m to cm
del_snowdepth = df_snowdepth.diff(periods=1,axis=1).ix[:,1]*100

# Move to fig dir
os.chdir(fig_dir)

# Plot on map
fig3, ax = make_map()
#kw = dict(marker='.', linestyle='-', alpha=0.25, color='darkgray', zorder=0)
#lines = ax.triplot(tri_info['X'], tri_info['Y'], **kw)
ax.set_title('Snowpack change.\n '+pd.to_datetime(time_stamp[0]).strftime('%Y-%m-%d %H:%M:%S MST')+' to '+pd.to_datetime(time_stamp[-1]).strftime('%Y-%m-%d %H:%M:%S MST'))
p1 = ax.tripcolor(tri_info['X'], tri_info['Y'], tri_info['triang'],facecolors=del_snowdepth,cmap='seismic_r',vmin=-50,vmax=50)
b1 = fig3.colorbar(p1)
b1.ax.set_ylabel('delta depth (cm)')
#sns.set_palette("BuGn_r")
fig3.savefig('48_snowdepth_change.png',bbox_inches='tight',dpi=300)



plt.show()
