import sys
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
sns.set_style('ticks')
sns.set_context("talk", font_scale=3, rc={"lines.linewidth": 2.5})
fig_res = 90

if len(sys.argv) == 1:
    sys.error('Missing arg of CHM run dir.')

chm_run_dir = str(sys.argv[-1])

# Dir to vtu files
main_dir  = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/SnowCast/CHM_Configs/')
vtu_dir   = os.path.join(main_dir,chm_run_dir,'meshes')
fig_dir   = os.path.join(main_dir,chm_run_dir,'figures')
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
ax.set_title('Snow depth.\n'+pd.to_datetime(time_stamp[0]).strftime('%Y-%m-%d %H:%M:%S MST')+'\n')
p1 = ax.tripcolor(tri_info['X'], tri_info['Y'], tri_info['triang'],facecolors=z*100,cmap='Blues',vmin=0, vmax=np.nanmax(z)*100*0.8) # Max at 80% of mean (helps show small depths) 
b1 = fig.colorbar(p1)
b1.ax.set_ylabel('Snow depth (cm)')
gl = ax.gridlines(draw_labels=True)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#sns.set_palette("BuGn_r")
fig.savefig('snowdepth_48h.png',bbox_inches='tight',dpi=fig_res)

# Input
var_name = 'swe'

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

# Get data
# Get mesh
cmesh = vfunc.get_mesh(cfile)
swe = vfunc.get_face_var(cmesh,var_name)
#df_swe = vfunc.get_multi_mesh_var_dask(cfile,var_name,time_stamp)

# Get swe, calc bulk density (kg m^2)
swe = swe/1000 # mm to m
bulk_density = swe/z*1000 # fraction to kg m^2
bulk_density[~np.isfinite(bulk_density)] = 0

# Move to fig dir
os.chdir(fig_dir)

# Plot on map
fig2, ax = make_map()
#kw = dict(marker='.', linestyle='-', alpha=0.25, color='darkgray', zorder=0)
#lines = ax.triplot(tri_info['X'], tri_info['Y'], **kw)
ax.set_title('Snow Density.\n '+pd.to_datetime(time_stamp[0]).strftime('%Y-%m-%d %H:%M:%S MST')+'\n')
p1 = ax.tripcolor(tri_info['X'], tri_info['Y'], tri_info['triang'],facecolors=bulk_density,cmap='Greens', vmin=50, vmax=450)
b1 = fig2.colorbar(p1)
b1.ax.set_ylabel('Bulk Density kg/m^3')
#sns.set_palette("BuGn_r")
fig2.savefig('density.png',bbox_inches='tight',dpi=fig_res)

## 48 accumued snow depth
# Input
var_name = 'snowdepthavg'

# Move to vtu dir
os.chdir(vtu_dir)
# Get list of files
vtu_files = np.sort(glob.glob(prefix+'*.vtu'))
# file to plot
cfile = [vtu_files[-2],vtu_files[-1]]

# Get time
time_stamp = vfunc.get_vtu_time(cfile,prefix)

# Adjust time zone
time_stamp = pd.to_datetime(time_stamp) + datetime.timedelta(hours=local_time_offset)

# Get data
df_snowdepth = vfunc.get_multi_mesh_var_dask(cfile,var_name,time_stamp)

# Take diff in 48 hrs m to cm
del_snowdepth = df_snowdepth.diff(periods=1,axis=1).ix[:,1]*100
abs_max = np.max([del_snowdepth.max(), -1*del_snowdepth.min()])

# Move to fig dir
os.chdir(fig_dir)

# Plot on map
fig3, ax = make_map()
#kw = dict(marker='.', linestyle='-', alpha=0.25, color='darkgray', zorder=0)
#lines = ax.triplot(tri_info['X'], tri_info['Y'], **kw)
ax.set_title('Snowpack change.\n '+pd.to_datetime(time_stamp[0]).strftime('%Y-%m-%d %H:%M:%S MST')+'\nto '+pd.to_datetime(time_stamp[-1]).strftime('%Y-%m-%d %H:%M:%S'))
p1 = ax.tripcolor(tri_info['X'], tri_info['Y'], tri_info['triang'],facecolors=del_snowdepth,cmap='seismic_r',vmin=-1*abs_max,vmax=abs_max)
b1 = fig3.colorbar(p1)
b1.ax.set_ylabel('delta depth (cm)')
#sns.set_palette("BuGn_r")
fig3.savefig('48_snowdepth_change.png',bbox_inches='tight',dpi=fig_res)



plt.show()
