import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import xarray as xr
import seaborn as sns
import sys
import imp


''' Plot spatial map at points, showing value for variable X'''
def plot_points_on_map():
