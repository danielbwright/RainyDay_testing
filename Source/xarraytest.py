#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 09:14:45 2023

@author: daniel
"""



def find_unique_elements(list1, list2):
    """
    Used to return only the elements of list1 that are not present in list2

    Parameters
    ----------
    list1 : target list of values to be reduced according to list2
    variables : list of values used to identify the values to keep in list1

    Returns
    -------
    list with only the values in list1 that were not present in list2

    """
    unique_elements_in_list1 = [x for x in list1 if x not in list2]
    unique_elements_in_list2 = [x for x in list2 if x not in list1]
    return unique_elements_in_list1



# invars is taken from 'variables' in the json file:
invars={"rainname": "precrate","latname":"latitude","longname":"longitude"}

# we don't want to drop these:
del invars['latname']   
del invars['longname']
keepvars=list(invars.values())

# open the "entire" netcdf file once in order to get the list of all variables:
infile=xr.open_dataset('/Users/daniel/Documents/DATA/StageIV/StageIV_InfilledCorr03degrees/FixedDate/StageIV.20020106.03degree.nc')

# this will only keep the variables that we need to read in. 
droplist=find_unique_elements(infile.keys(),keepvars) # droplist will be passed to the 'drop_variables=' in xr.open_dataset within the storm catalog creation loop in RainyDay
infile.close()

smallerfile=xr.open_dataset('/Users/daniel/Documents/DATA/StageIV/StageIV_InfilledCorr03degrees/FixedDate/StageIV.20020106.03degree.nc',drop_variables=droplist)



