#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 01:03:11 2023

@author: ashar
"""
import numpy as np
import xarray as xr
def readnetcdf(rfile,inbounds,variables=None):
    infile=xr.open_dataset(rfile)
    rain_name,lat_name,lon_name = variables.values()
    if np.any(inbounds!=False):
        latmin,latmax,longmin,longmax = inbounds[2],inbounds[3],inbounds[0],inbounds[1]
        outrain=infile[rain_name].sel(**{lat_name:slice(latmin,latmax)},\
                                                  **{lon_name:slice(longmin,longmax)})
        outlatitude=outrain[lat_name]
        outlongitude=outrain[lon_name]         
    else:
        outrain=infile[rain_name]
        outlatitude=outrain[lat_name]
        outlongitude=outrain[lat_name] 
    outtime=np.array(infile['time'],dtype='datetime64[m]')
    infile.close()
    return outrain,outtime,outlatitude,outlongitude

rfile = '/Volumes/TheCordex/Ben_Data/_1980/AORC.19800101.preciptemp.nc'
variables = {"rainarray":"RAINRATE", "latitude":"latitude", "longitude":"longitude"}
inarea = [-90.5, -88, 42.5, 44]
rain,timer,latr,lonr = readnetcdf(rfile, inarea, variables)
#%%
def readcatalog(rfile) :
    """
    Returns the properties of the storm including spatial range, storm center,
    storm depth, storm time by reading the already created storm catalogs.

    Parameters
    ----------
    rfile : string
        This takes in the path of the source file.

    Returns
    -------
    arrays
        This returns all the properties of a storm including storm rain array, storm time, storm depth, storm center and the extent of the transpositio domain

    """
    infile=xr.open_dataset(rfile)
    outrain=infile['rainrate']
    outlatitude=infile['latitude']
    outmask=infile['gridmask']
    domainmask=infile['domainmask']
    outtime=np.array(infile['time'],dtype='datetime64[m]')
    outlongitude=infile['longitude']
    outlocx=infile['xlocation']
    outlocy=infile['ylocation']
    outmax=infile['basinrainfall']

    try:
        timeresolution=np.int(infile.timeresolution)
        resexists=True
    except:
        resexists=False
    infile.close()
    
    if resexists:
        return outrain,outtime,outlatitude,outlongitude,outlocx,outlocy,outmax,outmask,domainmask,timeresolution
    else:
        return outrain,outtime,outlatitude,outlongitude,outlocx,outlocy,outmax,outmask,domainmask
