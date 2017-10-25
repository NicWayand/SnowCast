import sys
import os
import shutil 
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
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

def plot_variable(tri_var, var2plot, time_start, c_timestamp, fig_res, var_vmin, var_vmax):
    # Make map
    fig, ax = make_map()

    # Get time stamp as string
    str_timestamp = pd.to_datetime(c_timestamp).strftime('%Y-%m-%d %H:%M:%S MST')
    if time_start:
        str_time_start = pd.to_datetime(time_start).strftime('%Y-%m-%d %H:%M:%S')
    file_timestamp = pd.to_datetime(c_timestamp).strftime('%Y_%m_%d_%H:%M:%S_MST.png')
    
    # Make tri plot
    p1 = ax.tripcolor(tri_info['X'], tri_info['Y'], tri_info['triang'],
                      facecolors=tri_var * scale_factor[var2plot],
                      cmap=cmap_dict[var2plot],
                      vmin= var_vmin * scale_factor[var2plot],
                      vmax= var_vmax * scale_factor[var2plot])
    
    # Add title
    if time_start:
        ax.set_title(title_dict[var2plot] + '.\n' + str_time_start + ' to\n'+str_timestamp + '\n')
    else:
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
    if time_start:
        file_out = os.path.join(fig_dir, var2plot + '.png')
    else:
        file_out = os.path.join(fig_dir, var2plot, var2plot+'_'+file_timestamp)
    save_figure(fig, file_out, fig_res)
    return file_out

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

# Correct CHM units (NOT STANDARD METRIC....)
chm_units_fix = {'snowdepthavg':1, 'swe':1.0/1000} # to m

# Plotting dictionaries
ylabel_dict = {'snowdepthavg':'Snowdepth (cm)','swe':'SWE(mm)',
               'snowdepthavg_diff': 'Snowdepth change (cm)', 'swe_diff': 'SWE change (mm)'}
# From meteric (m) to display units
scale_factor = {'snowdepthavg':100,'swe':1000,
                'snowdepthavg_diff': 100, 'swe_diff': 1000,}
title_dict = {'snowdepthavg':'Snowdepth','swe':'SWE',
              'snowdepthavg_diff': 'Snowdepth change', 'swe_diff': 'SWE change'}
cmap_dict = {'snowdepthavg':mpl.colors.ListedColormap(sns.color_palette("Blues", 12)),
             'swe':mpl.colors.ListedColormap(sns.color_palette("Reds", 12)),
             'snowdepthavg_diff': mpl.colors.ListedColormap(sns.color_palette("RdBu", 12)),
             'swe_diff': mpl.colors.ListedColormap(sns.color_palette("RdBu", 12))}
var_min_delta = {'snowdepthavg':0.1,'swe':0.01} # Min max value for plotting change (m)


# Make single variable plots
var_names_2_plot = ['snowdepthavg','swe']
last_N = 2 # last output CHM vtu files to plot

# Make dirs
for var2plot in var_names_2_plot:
    if not os.path.isdir(os.path.join(fig_dir,var2plot)):
        os.mkdir(os.path.join(fig_dir,var2plot))

# Plot snapshots in time
for var2plot in var_names_2_plot:
    # Get data for all time stamps
    df_cvar = vfunc.get_multi_mesh_var_dask(ps.tail(last_N).values, var2plot, ps.tail(last_N).index)
    df_cvar = df_cvar * chm_units_fix[var2plot] # Units to Metric standard (m)
    df_cvar_max = df_cvar.max().max()*0.8 # Set max at 80%
    df_cvar_min = 0
    for ct in df_cvar.index:
        print('printing ',ct)
        # Plot and Save
        file_out = plot_variable(df_cvar.loc[ct].values, var2plot, [], ct, fig_res, df_cvar_min, df_cvar_max)
        print('fig saved')
        
        # If last time
        if ct==df_cvar.index[-1]:
            print('Copied last vtu to '+var2plot+'_most_recent.png')
            shutil.copyfile(file_out, os.path.join(fig_dir,var2plot+'_most_recent.png'))


# Plot change from NOW* to most future forecast
# *NOW is defined as the CHM output time stamp most nearest to the time plot sripts run
time_now = datetime.datetime.now() - datetime.timedelta(hours=1)# Chinook is on CST (-6) need to subtract 1 hour to got MST -7
i = np.argmin(np.abs(ps.index.to_pydatetime() - time_now))
time_now_near = ps.index[i]
time_fut = ps.index[-1]

for var2plot in var_names_2_plot:
    # Get data for all time stamps
    df_cvar = vfunc.get_multi_mesh_var_dask(ps.loc[[time_now_near, time_fut]].values,
                                            var2plot, ps.loc[[time_now_near, time_fut]].index)
    df_cvar = df_cvar * chm_units_fix[var2plot] # Units to Metric standard (m)
    # Take difference
    df_dif = df_cvar.diff(axis=0).iloc[1] # diff over time, returns missing on first row, so take second
    df_trim = df_dif[np.abs(df_dif-df_dif.mean())<=(3*df_dif.std())] # Remove outliers (sometimes from avalnching)
    df_cvar_max = np.max([np.abs(df_trim.min()), df_trim.max()]) # Use larger of neg or pos values
    df_cvar_max = np.max([df_cvar_max, var_min_delta[var2plot]]) # Min max delta value
    df_cvar_min = -1*df_cvar_max
    print(df_cvar_min, df_cvar_max)
    print('printing ',time_fut,' minus ',time_now_near)
    # Plot and Save
    file_out = plot_variable(df_dif.values, var2plot+'_diff', time_now_near, time_fut, fig_res, df_cvar_min, df_cvar_max)
    print('fig saved')

import sys
sys.exit()
### OLD FIGURE CODE BELOW

cfile = vtu_files[-1]

# Get time
time_stamp = vfunc.get_vtu_time([cfile],prefix)

# Adjust time zone
time_stamp = pd.to_datetime(time_stamp) + datetime.timedelta(hours=local_time_offset)

# Get data
# Get mesh
cmesh = vfunc.get_mesh(cfile)

z = vfunc.get_face_var(cmesh,'snowdepthavg')

# Move to fig dir
os.chdir(fig_dir)

# Plot on map
fig, ax = make_map()
#kw = dict(marker='.', linestyle='-', alpha=0.25, color='darkgray', zorder=0)
#lines = ax.triplot(tri_info['X'], tri_info['Y'], **kw)
ax.set_title('Snow depth.\n'+pd.to_datetime(time_stamp[0]).strftime('%Y-%m-%d %H:%M:%S MST')+'\n')
p1 = ax.tripcolor(tri_info['X'], tri_info['Y'], tri_info['triang'],facecolors=z*100,cmap=cmap_dict[var2plot],vmin=0, vmax=np.nanmax(z)*100*0.8) # Max at 80% of mean (helps show small depths)
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
