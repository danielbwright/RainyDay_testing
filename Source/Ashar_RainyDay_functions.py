import os
import sys
import numpy as np
import xarray as xr
from numpy import unravel_index
import scipy as sp
import glob
import math
from datetime import datetime, date, time, timedelta
import time
from copy import deepcopy
import fiona

import cartopy
from matplotlib.patches import Polygon
from scipy import stats
from netCDF4 import Dataset, num2date, date2num
import rasterio
from rasterio.transform import from_origin
from rasterio.shutil import delete
from rasterio.mask import mask
from rasterio.io import MemoryFile
import pandas as pd
from numba import prange, jit
import dask.array as da
import json
import re
import pyproj
from shapely.ops import transform

import shapely
from shapely.ops import unary_union
from shapely.geometry import shape

import geopandas as gp

from scipy.stats import norm
from scipy.stats import lognorm

# plotting stuff, really only needed for diagnostic plots
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm


###Pass the dataset from xarray and output will be dataset with excluded years and month.
def createfilelist(ds, includeyears, excludemonths):
    if isinstance(includeyears, (bool)) == False:
        ds = ds.loc[{"time": ds.time.dt.year.isin(includeyears)}]

    if isinstance(excludemonths, (bool)) == False:
        ds = ds.loc[{"time": ~ds.time.dt.month.isin(excludemonths)}]

    nyears = len(np.unique(ds.time.dt.year))
    return ds, nyears


def readnetcdf(ds, inbounds=False, lassiterfile=False):
    # infile =ds
    if lassiterfile == False:
        if 'rainrate' in ds.variables.keys():

            print("Reading an older-style NetCDF file")
            oldfile = True
        else:
            oldfile = False
            precrate = ds.precrate

        if np.any(inbounds != False):
            latitude = ds.latitude
            longitude = ds.longitude
            time = ds.time
            rain = ds.rainrate
            if oldfile:
                outrain = rain[:, inbounds[3]:inbounds[2] + 1, inbounds[0]:inbounds[1] + 1]
                outlatitude = latitude[inbounds[3]:inbounds[2] + 1]
            else:
                outrain = precrate[:, ::-1, :][:, inbounds[3]:inbounds[2] + 1, inbounds[0]:inbounds[1] + 1]
                outlatitude = latitude[::-1][inbounds[3]:inbounds[2] + 1]
            outlongitude = longitude[inbounds[0]:inbounds[1] + 1]
        else:
            day_ds = ds.isel(time=slice(0, 24))
            latitude = day_ds.latitude
            longitude = day_ds.longitude
            time = day_ds.time
            rain = day_ds.rainrate
            if oldfile:
                outrain = rain
                outlatitude = latitude
            else:
                outrain = precrate[:, ::-1, :]
                outlatitude = latitude[::-1]
            outlongitude = longitude
        outtime = time.values.astype('datetime64[m]')
    else:  # lassiter time!
        print("Lassiter  or FitzGerald style!")
        # for subhourly lassiter files:
        if np.any(inbounds != False):
            outrain = rain[:, ::-1, :][:, inbounds[3]:inbounds[2] + 1, inbounds[0]:inbounds[1] + 1]
            outlatitude = latitude[inbounds[3]:inbounds[2] + 1][::-1]
            outlongitude = longitude[:] - 360.
        else:
            outrain = rain[:][:, ::-1, :]
            outlatitude = latitude[:][::-1]
            outlongitude = longitude[:] - 360.

            # for hourly lassiter files:
        # tempdate=rfile.strip('.nc').split('/')[-1]
        # startdate=np.datetime64(tempdate[0:4]+'-'+tempdate[4:6]+'-'+tempdate[6:8]+'T00:00')
        # outtime=startdate+np.array(infile.variables['time'][:],dtype='timedelta64[s]')

        # for hourly FitzGerald files:
        #### Ashar: Might not work for now.. ####
        tempdate = rfile.strip('.nc').split('/')[-1][-8:]
        startdate = np.datetime64(tempdate[0:4] + '-' + tempdate[4:6] + '-' + tempdate[6:8] + 'T00:00')
        outtime = startdate + np.array(time[:], dtype='timedelta64[m]')

        # if np.any(inbounds!=False):
        #   outrain=np.array(infile.variables['precrate'][::-1][:,inbounds[3]:inbounds[2]+1,inbounds[0]:inbounds[1]+1])
        #    outlatitude=np.array(infile.variables['latitude'][::-1][inbounds[3]:inbounds[2]+1])
        #    outlongitude=np.array(infile.variables['longitude'][inbounds[3]:inbounds[2]+1])
        # else:
        #   outrain=np.array(infile.variables['precrate'][:,::-1,:])
        #    outlatitude=np.array(infile.variables['latitude'][::-1])
        #    outlongitude=np.array(infile.variables['longitude'][:])

    return outrain, outtime, outlatitude[1], outlongitude[1]


def rainprop_setup(infile, catalog=False, lassiterfile=False):
    if catalog:
        inrain, intime, inlatitude, inlongitude, catx, caty, catmax, _, domainmask = readcatalog(infile)
    else:
        inrain, intime, inlatitude, inlongitude = readnetcdf(infile, lassiterfile=lassiterfile)

    if len(inlatitude.shape) > 1 or len(inlongitude.shape) > 1:
        sys.exit("RainyDay isn't set up for netcdf files that aren't on a regular lat/lon grid!")
        # inlatitude=inlatitude[:,0]          # perhaps would be safer to have an error here...
        # inlongitude=inlongitude[0,:]        # perhaps would be safer to have an error here...
    yres = (inlatitude[0:-1] - inlatitude[1:]).mean()
    xres = (inlongitude[1:] - inlongitude[0:-1]).mean()
    if np.isclose(xres, yres) == False:
        sys.exit("Rainfall grid isn't square. RainyDay cannot support that.")

    unqtimes = np.unique(intime)
    if len(unqtimes) > 1:
        tdiff = unqtimes[1:] - unqtimes[0:-1]
        tempres = np.min(unqtimes[1:] - unqtimes[0:-1])  # temporal resolution
        if np.any(np.not_equal(tdiff, tempres)):
            sys.exit("Uneven time steps. RainyDay can't handle that.")
    else:
        # this is to catch daily data where you can't calculate a time resolution
        tempres = np.float32(1440)
        tempres = tempres.astype(
            'timedelta64[m]')  # temporal resolution in minutes-haven't checked to make sure this works right

    if len(intime) * np.float32(tempres) != 1440. and catalog == False:
        sys.exit("RainyDay requires daily input files, but has detected something different.")
    tempres = int(np.float32(tempres))

    nodata = inrain.groupby(inrain.where(inrain < 0.)).count()
    if len(nodata) > 1:
        sys.exit("More than one missing value flag.")
    elif len(nodata) == 0 and catalog == False:
        print(
            "Warning: Missing data flag is ambiguous. RainyDay will probably handle this ok, especially if there is not missing data.")
        nodata == -999.
    elif catalog:
        nodata = -999.
    else:
        nodata = nodata[0]

    if catalog:
        return [xres, yres], [len(inlatitude), len(inlongitude)], [np.min(inlongitude), np.max(inlongitude),
                                                                   np.min(inlatitude), np.max(
                inlatitude)], tempres, nodata, inrain, intime, inlatitude, inlongitude, catx, caty, catmax, domainmask
    else:
        return [xres, yres], [len(inlatitude), len(inlongitude)], [inlongitude.min(), inlongitude.max() + xres,
                                                                   inlatitude.min() - yres,
                                                                   inlatitude.max()], tempres, nodata


def findsubbox(inarea, rainprop):           ### Taking a lot of time with dask.
    outind = da.empty([4], dtype='int')
    # outextent=np.empty([4])
    outdim = da.empty([2])
    res = rainprop.spatialres[0].values
    # inbox=[inarea[0]+res/2.,inarea[1]+res/2.,inarea[2]-res/2.,inarea[3]-res/2.]

    # DBW: updated the stuff below on 9/13/22 to accomodate cell-centers rather than upper-left corners. Older files will probably have some issues here.....
    rangex = da.arange(rainprop.bndbox[0], rainprop.bndbox[1] - res / 1000, res)
    rangey = da.arange(rainprop.bndbox[3], rainprop.bndbox[2] + res / 1000, -res)

    if rangex.shape[0] < rainprop.dimensions[1]:
        rangex = np.append(rangex, rangex[-1])
    if rangey.shape[0] < rainprop.dimensions[0]:
        rangey = np.append(rangey, rangey[-1])
    if rangex.shape[0] > rainprop.dimensions[1]:
        rangex = rangex[0:-1]
    if rangey.shape[0] > rainprop.dimensions[0]:
        rangey = rangey[0:-1]

    outextent = inarea

    # "SNAP" output extent to grid
    outind[0] = da.abs((rangex - res / 2.) - outextent[0]).argmin()
    outind[1] = da.abs((rangex + res / 2.) - outextent[1]).argmin() - 1
    outind[2] = da.abs((rangey - res / 2.) - outextent[2]).argmin() - 1
    outind[3] = da.abs((rangey + res / 2.) - outextent[3]).argmin()
    outextent[0],outextent[1] = rangex[outind[0]], rangex[outind[1] + 1]
    outextent[2], outextent[3] = rangey[outind[2] + 1], rangey[outind[3]]

    outdim[1]=da.shape(da.arange(outind[0],outind[1]+1))[0]
    outdim[0]=da.shape(da.arange(outind[3],outind[2]+1))[0]
    outdim = np.array(outdim, dtype='int32')
    return outextent, outind, outdim


@jit(parallel=True)
def catalogNumba(temparray, trimmask, xlen, ylen, maskheight, maskwidth, rainsum):
    rainsum[:] = 0.
    for i in range(0, (ylen) * (xlen)):
        y = i // xlen
        x = i - y * xlen
        # print x,y
        rainsum[y, x] = da.nanmean(
            temparray[(y):(y + maskheight), (x):(x + maskwidth)] * trimmask)  ### Take the mean of rainfall at a region.

    # wheremax=np.argmax(rainsum)
    rainsum  = da.where(rainsum < 5, 0, rainsum)
    rmax     = da.nanmax(rainsum)
    wheremax = da.where(da.equal(rainsum, rmax))
    return rmax, wheremax[0].compute()[0], wheremax[1].compute()[0]


def writecatalog(scenarioname, catrain, catmax, catx, caty, cattime, latrange, lonrange, catalogname, nstorms, gridmask,
                 parameterfile, dmask, timeresolution=False):

    with open(parameterfile,'r') as f:
        cardinfo = json.loads(f.read())
    # Variable Attributes (time since 1970-01-01 00:00:00.0 in numpys)
    latitudes_units,longitudes_units = 'degrees_north', 'degrees_east'
    rainrate_units,basinrainfall_units = 'mm hr^-1', 'mm'
    times_units, times_calendar = 'minutes since 1970-01-01 00:00.0' , 'gregorian'
    
    # Variable Names
    times_name = 'time'                     ## change here
    latitudes_name,longitudes_name  = 'latitude', 'longitude'
    rainrate_name, basinrainfall_name = 'precipitation rate', 'storm total basin averaged precipitation'
    xlocation_name, ylocation_name = 'x index of storm', 'y index of storm'
    gmask_name, domainmask_name = 'mask for Aw (control volume)', 'mask for transposition domain'

    history, missing = 'Created ' + str(datetime.now()), '-9999.'
    source = 'RainyDay Storm Catalog for scenario ' + scenarioname + '. See description for JSON file contents.'
    data_vars = dict(precrate = (("time","longitude","latitude"),catrain[:, :, ::-1, :],{'units': rainrate_units, 'long_name': rainrate_name}),
                     basinrainfall = ((),catmax,{'units': basinrainfall_units, 'long_name': basinrainfall_name}),
                     xlocation = ((),catx,{'units': 'dimensionless', 'long_name': xlocation_name}),
                     ylocation = ((),catx,{'units': 'dimensionless', 'long_name': xlocation_name}),
                     gridmask= (("latitude", "longitude"), gridmask[::-1, :],{'units': 'dimensionless', 'long_name': gmask_name}),
                     domainmask = (("latitude", "longitude"),dmask[::-1, :],{'units': 'dimensionless', 'long_name': domainmask_name}))
    coords = dict(time = (("time"),cattime,{'units':times_units,'long_name':times_name}),
                  longitude = (("longitude"), lonrange, {'units': longitudes_units, 'long_name': longitudes_name}),
                  latitude =  (("latitude"), latrange[::-1], {'units': latitudes_units, 'long_name': latitudes_name}))


    attrs  = dict(history =history, source =  source, missing = missing, description = cardinfo,  calendar = times_calendar)
    catalog = xr.Dataset(data_vars = data_vars, coords = coords, attrs = attrs)
    catalog.to_netcdf(catalogname)




 


































