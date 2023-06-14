#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 01:03:11 2023

@author: ashar
"""
import numpy as np
import xarray as xr
def readnetcdf(rfile,inbounds=False,lassiterfile=False):
    infile=xr.open_dataset(rfile)
    search_string = 'rainrate'     #### Used these lines to counter small case letters.
    variables = list(infile.variables.keys())
    index = [x.lower() for x in variables].index(search_string.lower())
    if np.any(inbounds!=False):
        latmin,latmax,longmin,longmax = inbounds[2],inbounds[3],inbounds[0],inbounds[1]
        outrain=infile[variables[index]].sel(latitude =slice(latmin,latmax),\
                                                  longitude=slice(longmin,longmax))
        outlatitude=outrain['latitude']
        outlongitude=outrain['longitude']        
    else:
        outrain=infile[variables[index]]
        outlatitude=outrain['latitude']
        outlongitude=outrain['longitude'] 
    outtime=np.array(infile['time'],dtype='datetime64[m]')
    infile.close()
    return outrain,outtime,outlatitude,outlongitude

rfile = '/Volumes/TheCordex/Ben_Data/_1980/AORC.19800101.preciptemp.nc'
inarea = [-90.5, -88, 42.5, 44]
rain,timer,latr,lonr = readnetcdf(rfile, inarea)