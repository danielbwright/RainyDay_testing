#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 20:56:59 2023

@author: ashar
"""

import zipfile
import os
import xarray as xr
import numpy as np
import sys

# Specify the path to the zip file
storm_path = '/Volumes/TheCordex/Madison_Data/'
storm_yr = np.array([])
storms = int(sys.argv[1])
for i in range(storms):
    storm = xr.open_dataset(storm_path + 'Storm'+str(i+1)+'.nc')
    start_time = int(storm.time[0].dt.year)
    storm_yr = np.append(storm_yr, start_time)

unique, ncounts = np.unique(storm_yr, return_counts=True)
count_storms = dict(zip(unique, ncounts))

print(count_storms)

        