#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 20:47:06 2023

@author: ashar
"""

import numpy as np
from datetime import datetime
import glob
import pandas as pd
import os
# Create a NumPy datetime variable
date_num =np.datetime64('1980-01-21T00:00:00')
print(date_num + 7200)
# datetime_var =  datetime.strptime('2008-01-01T12:34:56', "%Y-%m-%dT%H:%M:%S")
# print(date_num.date())
#%%
# Print the datetime variable
datafolder = "/Volumes/TheCordex/CONUS404/1980/"
current = date_num.date()
file_pattern = datafolder + current.strftime("*%Y-%m-%d.nc")
file_paths = glob.glob(file_pattern)

#%%
nstorms = 400
for i in range(nstorms):
    start_time = cattime[i,0]
    end_time   = cattime[i,-1]
    current_datetime = start_time
    dataset_date = None
    catrain = np.zeros((int(catduration*60/rainprop.timeres),rainprop.subdimensions[0]\
                        ,rainprop.subdimensions[1]),dtype='float32')
    k  = 0 
    infile = None
    while current_datetime <= end_time:
        current_date = np.datetime_as_string(current_datetime, unit='D')
        if current_date != dataset_date:
            datset_date = current_date
            # This loop searches for the right file to open from the filelist
            for file in flist:
                match = re.search(r'\d{4}(?:\d{2})?(?:\d{2}|\-\d{2}\-\d{2}|\/\d{2}/\d{2})',\
                                  os.path.basename(file))
                if current_date.replace("-","") == match.group().replace("-","").replace("/",""):
                    infile = file
                    break
            stm_rain,stm_time,_,_ = RainyDay.readnetcdf(infile,inbounds=rainprop.subind,\
                                                        lassiterfile=islassiter)
        cind = np.where(stm_time == current_datetime)[0][0]
        catrain[k,:] = stm_rain[cind,:]
        current_datetime += rainprop.timeres * 60
        k += 1
    RainyDay.writecatalog_ash(scenarioname,catrain,\
                              catmax[i],\
                                  catx[i],caty[i],\
                                      cattime[i,:],catrain,latrange,lonrange,\
                                          storm_name,catmask,parameterfile_json,domainmask,\
                                              timeresolution=rainprop.timeres)
    
#%%
stime = np.datetime64('1980-01-21T00:00:00', 'm'); etime = np.datetime64('1980-01-21T12:00:00')
timestep = np.timedelta64(60,'m')
time_array = np.arange(stime,etime,timestep)
print(np.where(time_array == stime)[0][0])
#%%
for i in range(nstorms):
    start_time = cattime[i,0]
    end_time   = cattime[i,-1]
    current_datetime = start_time
    dataset_date = None
    catrain = np.zeros((int(catduration*60/rainprop.timeres),rainprop.subdimensions[0]\
                        ,rainprop.subdimensions[1]),dtype='float32')
    k  = 0 
    infile = None
    while current_datetime <= end_time:
        current_date = np.datetime_as_string(current_datetime, unit='D')
        if current_date != dataset_date:
            datset_date = current_date
            # This loop searches for the right file to open from the filelist
            for file in flist:
                match = re.search(r'\d{4}(?:\d{2})?(?:\d{2}|\-\d{2}\-\d{2}|\/\d{2}/\d{2})',\
                                  os.path.basename(file))
                if current_date.replace("-","") == match.group().replace("-","").replace("/",""):
                    infile = file
                    break
            stm_rain,stm_time,_,_ = RainyDay.readnetcdf(infile,inbounds=rainprop.subind,\
                                                        lassiterfile=islassiter)
        cind = np.where(stm_time == current_datetime)[0][0]
        catrain[k,:] = stm_rain[cind,:]
        current_datetime += rainprop.timeres * 60
        k += 1
    storm_name = "Storm" + str(i+1) +".nc"
    print("Writing Storm "+ str(i+1) + "out of " + str(nstorms) )
    RainyDay.writecatalog_ash(scenarioname,catrain,\
                              catmax[i],\
                                  catx[i],caty[i],\
                                      cattime[i,:],latrange,lonrange,\
                                          storm_name,catmask,parameterfile_json,domainmask,\
                                              timeresolution=rainprop.timeres)
#%%
stprm1 = xr.open_dataset("/Volumes/TheCordex/Madison_Data/Storm1.nc")
stprm2 = xr.open_dataset("/Volumes/TheCordex/Madison_Data/Storm199.nc")
stprm3 = xr.open_dataset("/Volumes/TheCordex/Madison_Data/Storm200.nc")