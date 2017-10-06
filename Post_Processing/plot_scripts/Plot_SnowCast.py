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
import seaborn as sns

# General plotting settings
sns.set_style('ticks')
sns.set_context("talk", font_scale=3, rc={"lines.linewidth": 2.5})
fig_res = 90 # dpi

if len(sys.argv) == 1:
    sys.error('Missing arg of CHM run dir.')

chm_run_dir = str(sys.argv[-1])

# Dir to vtu files
main_dir  = os.path.normpath(r'/home/nwayand/snow_models/output_CHM/SnowCast/CHM_Configs/')
vtu_dir   = os.path.join(main_dir,chm_run_dir,'meshes')
fig_dir   = os.path.join(main_dir,chm_run_dir,'figures')
prefix = 'SC'

local_time_offset = -7 # Offset to local standard time (i.e. -7 for MST)

# Make fig dir
if not os.path.isdir(fig_dir):
    os.mkdir(fig_dir)

# General plotting functions

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

def save_figure(f,file_out,fig_res):
    f.savefig(file_out,bbox_inches='tight',dpi=fig_res)    

def plot_variable(tri_var, var2plot, c_timestamp, fig_res, var_vmax):
    # Make map
    fig, ax = make_map()

    # Get time stamp as string
    str_timestamp = pd.to_datetime(c_timestamp).strftime('%Y-%m-%d %H:%M:%S MST')
    file_timestamp = pd.to_datetime(c_timestamp).strftime('%Y_%m_%d_%H:%M:%S_MST.png')
    
    # Make tri plot
    p1 = ax.tripcolor(tri_info['X'], tri_info['Y'], tri_info['triang'], facecolors=tri_var * scale_factor[var2plot],cmap='Blues',vmin=np.nanmin(tri_var), vmax= var_vmax * scale_factor[var2plot])
    
    # Add title
    ax.set_title(title_dict[var2plot]+'.\n'+str_timestamp+'\n')

    # Add colorbar
    b1 = fig.colorbar(p1)
    b1.ax.set_ylabel(ylabel_dict[var2plot])
    
    # Add grid lines
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Save figure
    file_out = os.path.join(fig_dir, var2plot, var2plot+'_'+file_timestamp) 
    save_figure(fig, file_out, fig_res)

# Get General info about vtus

# Move to vtu dir
os.chdir(vtu_dir)
# Get list of files
vtu_files = np.sort(glob.glob(prefix+'*.vtu'))
# Get time stamps
time_stamps_all_UTC = vfunc.get_vtu_time(vtu_files, prefix)
# Adjust time zone
time_stamps_all_LST = time_stamps_all_UTC + pd.DateOffset(hours=local_time_offset)
# Create series of time stamps and vtu files
ps = pd.Series(vtu_files, index=time_stamps_all_LST)
# Get mesh
cmesh = vfunc.get_mesh(vtu_files[-1])
# Get tri info
tri_info = vfunc.get_triangle(cmesh)

# Plotting dictionaries
ylabel_dict = {'snowdepthavg':'Snowdepth (cm)','swe':'SWE (mm)'}
scale_factor = {'snowdepthavg':100,'swe':1.0}
title_dict = {'snowdepthavg':'Snowdepth','swe':'SWE'}

# Make single variable plots
var_names_2_plot = ['snowdepthavg','swe']
last_N = 10 # last output CHM vtu files to plot

# Make dirs
for var2plot in var_names_2_plot:
    if not os.path.isdir(os.path.join(fig_dir,var2plot)):
        os.mkdir(os.path.join(fig_dir,var2plot))

for var2plot in var_names_2_plot:
    # Get data for all time stamps
    df_cvar = vfunc.get_multi_mesh_var_dask(ps.tail(last_N).values, var2plot, ps.tail(last_N).index)
    df_cvar_max = df_cvar.max().max()*0.8 # Set max at 80%
    for ct in df_cvar.index:
        # Plot and Save
        plot_variable(df_cvar.loc[ct].values, var2plot, ct, fig_res, df_cvar_max)
        print('fig saved')

# Copy last snowdepth and swe to orig files

import sys
sys.exit()

z = vfunc.get_face_var(cmesh,var_name)

# Move to fig dir
os.chdir(fig_dir)

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
# Then get times for ALL files
time_stamp_all = vfunc.get_vtu_time(vtu_files,prefix)
# last_time = time_stamp_all[-1]
t_dt = (time_stamp_all[1]-time_stamp_all[0])
time_dt_h = t_dt.days * 24 +  t_dt.seconds/3600
day_before_index = -(1+ int(24/time_dt_h))

# file to plot
cfile = [vtu_files[day_before_index],vtu_files[-1]]

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
