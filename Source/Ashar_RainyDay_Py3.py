import os
import sys
import numpy as np
import json
import scipy as sp
import shapefile
import math
from datetime import datetime, date, time, timedelta
import time
from copy import deepcopy
import dask.array as da
from scipy import ndimage, stats
import pandas as pd

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.io.shapereader import Reader
import glob
import geopandas as gp
import xarray as xr
from cartopy.feature import ShapelyFeature
import cartopy.mpl.ticker as cticker
#from numba import njit, prange
numbacheck=True

# plotting stuff, really only needed for diagnostic plots
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

# import RainyDay functions
import RainyDay_utilities_Py3.Ashar_RainyDay_functions as RainyDay

from numba.types import int32

import warnings
warnings.filterwarnings("ignore")

import tracemalloc
tracemalloc.start()
#==============================================================================
# RAINFALL CLASS
# THIS CONTAINS INFORMATION ABOUT THE SPECIFIED INPUT RAINFALL DATASET
#==============================================================================
class GriddedRainProperties(object):
    def __init__(self,dataset,bndbox,subind,subextent,dimensions,subdimensions,spatialres,timeres,timeunits,spatialunits,rainunits,nodata,notes):
        self.dataset=dataset                #  INPUT RAINFALL DATA SOURCE (TMPA, STAGE IV, ETC.)
        self.bndbox=bndbox                  # COORDINATES OF FULL DATASET BOUNDING BOX (THIS SHOULD BE THE IN THE ORDER [WEST LON,EAST LON,SOUTH LAT,NOTH LAT])
        self.subind=subind                  # MATRIX INDICIES OF BOUNDING BOX (XMIN,XMAX,YMIN,YMAX)
        self.subextent=subextent            # COORDINATES OF USER-DEFINED BOUNDING BOX (NOTE: THIS WILL BE AUTOMATICALLY CALCULATED)
        self.dimensions=dimensions          # WHAT ARE THE SPATIAL DIMENSIONS OF THE INPUT DATASET (NOTE: EVEN IF THE INPUT DATA YOU ARE USING IS FOR A SUBDOMAIN, THIS SHOULD BE THE DIMENSIONS OF THE FULL DATASET)
        self.subdimensions=subdimensions    # WHAT ARE THE DIMENSIONS OF THE SUBDOMAIN (IN THIS SCRIPT, THIS WILL BE RESET TO THE ACTUAL DIMENSION OF THE INPUT DATA, IF NEEDED)
        self.spatialres=spatialres          # WHAT IS THE SPATIAL RESOLUTION OF THE INPUT DATA [dx,dy]?  CURRENTLY THIS SCRIPT WILL ONLY HANDLE RECTANGULAR GRIDS IN DEGREES
        self.timeunits=timeunits            # TEMPORAL UNITS.  CURRENTLY MUST BE MINUTES
        self.spatialunits=spatialunits      # SPATIAL UNITS (CURRENTLY MUST BE DEGREES) [Xres,Yres]
        self.rainunits=rainunits            # RAINFALL UNITS (CURRENTLY MUST BE MM/HR)
        self.nodata=nodata                  # MISSING DATA FLAG
        self.notes=notes                    # ANY SPECIAL NOTES?


#==============================================================================
# RAINFALL INFO
# NOTE: "BOUNDING BOXES"-bndbox, subbox, CONSIDER THE COORDINATES OF THE EDGE OF THE BOUNDING BOX
#==============================================================================

emptyprop=GriddedRainProperties('emptyprop',
                            [-999.,-999.-999.-999.],
                            [999, 999, 999, 999],
                            [-999.,-999.-999.-999.],
                            [999, 999],
                            [999, 999],
                            [99.,99.],
                            99.,
                            "minutes",
                            "degrees",
                            "mm/hr",
                            -9999.,
                            "none")


################################################################################
# "MAIN"
################################################################################

print('''Welcome to RainyDay, a framework for coupling remote sensing precipitation
        fields with Stochastic Storm Transposition for assessment of rainfall-driven hazards.
        Copyright (C) 2017  Daniel Benjamin Wright (danielb.wright@gmail.com)
        Distributed under the MIT Open-Source License: https://opensource.org/licenses/MIT
    
    
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
''')


#==============================================================================
# USER DEFINED INFO
#==============================================================================

start = time.time()
parameterfile='ttt'
try:
    parameterfile=str(sys.argv[1])
    print("reading in the parameter file...")
    #### Cardinfo takes in the  '.sst' file parameters
    with open(parameterfile, 'r') as read_file:
        cardinfo = json.loads(read_file.read())
    # cardinfo=np.loadtxt(parameterfile, comments="#",dtype="str", unpack=False)
    #parameterfile='/Users/daniel/Google_Drive/Presentations/MyPresentations/WisconsinRainfallProject/Wisconsin_24hrIDF_Madison.sst'
except :
    print("You either didn't specify a parameter file, or it doesn't exist on the source path given.")

#### Changed this part of the code
# if os.path.isfile(parameterfile)==False:
#     sys.exit("You either didn't specify a parameter file, or it doesn't exist.")
# else:
#     print("reading in the parameter file...")
#     #### Cardinfo takes in the  '.sst' file parameters
#     cardinfo=np.loadtxt(parameterfile, comments="#",dtype="str", unpack=False)
#
#
# #==============================================================================
# # USER DEFINED VARIABLES
# #==============================================================================
#
try:
    wd=cardinfo["MAINPATH"]
    #wd=np.str(sys.argv[2])
    try:
        scenarioname=cardinfo["SCENARIONAME"]
    except ValueError:
        print("You didn't specify SCENARIONAME, defaulting to 'RainyDayAnalysis'!")
        scenarioname='RainyDayAnalysis'
    fullpath=wd+'/'+scenarioname
    if os.path.isdir(fullpath)==False:
        #os.system('mkdir -p -v %s' %(fullpath))
        os.mkdir(fullpath)
    os.chdir(wd)
#### Changed this part of the code
except ImportError:
    sys.exit("You did't specify MAINPATH, which is a required field!")
#
#
# # PROPERTIES RELATED TO SST PROCEDURE:
try:
    catalogname=cardinfo["CATALOGNAME"]
except Exception:
    catalogname='SSTstormcatalog.nc'
if catalogname.find('/')==-1:
    catalogname=wd+'/'+catalogname     ### This doesn't looks good to me. We can avaoid this "if" and jsut add the path.

if len(catalogname.split('.nc'))==1:
    sys.exit("You didn't specify the '.nc' extension on the storm catalog name!")
#
try:
    CreateCatalog=cardinfo["CREATECATALOG"]
except Exception:
    sys.exit("You didn't specify CREATECATALOG, which is required!")

if CreateCatalog.lower()=='true':
    CreateCatalog=True
else:
    CreateCatalog=False
    if os.path.isfile(catalogname)==False:
        sys.exit("You need to create a storm catalog first.")
    else:
        print("Reading an existing storm catalog!")
        try:
            catrain,cattime,latrange,lonrange,catx,caty,catmax,_,domainmask,timeres=RainyDay.readcatalog(catalogname)
        except ValueError:
            catrain,cattime,latrange,lonrange,catx,caty,catmax,_,domainmask=RainyDay.readcatalog(catalogname)
            timeres=RainyDay.readtimeresolution(catalogname)
        #if catrain.shape[1]==1:
        #    timeres=RainyDay.readtimeresolution(catalogname)
        yres=np.abs(np.mean(latrange[1:]-latrange[0:-1]))
        xres=np.abs(np.mean(lonrange[1:]-lonrange[0:-1]))
        catarea=[lonrange[0],lonrange[-1]+xres,latrange[-1]-yres,latrange[0]]
        if np.isclose(xres,yres)==False:
            sys.exit('Sadly, RainyDay currently only supports equal x and y resolutions.')
        else:
            res=np.min([yres,xres])
#
#
try:
    nstorms=cardinfo["NSTORMS"]
    defaultstorms=False
    print(str(nstorms)+" storms will be identified for storm catalog!")
except Exception:
    #nstorms=100
    print("you didn't specify NSTORMS, defaulting to 20 per year, or whatever is in the catalog!")
    defaultstorms=True
#
try:
    nsimulations=cardinfo["NYEARS"]
    print(str(nsimulations)+" years of synthetic storms will be generated, for a max recurrence interval of "+str(nsimulations)+" years!")
except Exception:
    nsimulations=100
    print("you didn't specify NYEARS, defaulting to 100, for a max recurrence interval of 100 years!")
#
try:
    nrealizations=cardinfo["NREALIZATIONS"]
    print(str(nrealizations)+" realizations will be performed!")
except Exception:
    nrealizations=1
#
#
try:
    duration=cardinfo["DURATION"]
    if duration<=0.:
        sys.exit("Duration is zero or negative!")
except ImportError:
    sys.exit("You didn't specify 'DURATION', which is a required field!")
#
#
# # this following bit is only needed in the very specific (and generally not recommended) case when the desired analysis/catalog duration is equal to the temporal resolution of the dataset
# # try:
# #     temptimeres=np.int(cardinfo[cardinfo[:,0]=="TIMERESOLUTION",1][0])
# #     print("A resolution of "+str(temptimeres)+" minutes has been provided. Be careful with this, because if it is improperly specified, this will cause errors. Note that TIMERESOLUTION is not needed unless the duration of each storm is to be exactly equal to the temporal resolution of the input data or catalog. In other words, make sure you know what you're doing!")
# # except Exception:
# #     temptimeres=False
#
#
#
try:
    samplingtype=cardinfo["RESAMPLING"]
    sampling_options = ['poisson', 'empirical', 'negbinom']
    if samplingtype.lower() not in sampling_options:
        sys.exit("unrecognized storm count resampling type, options are: 'poisson', 'empirical', or 'negbinom'")
except Exception:
    samplingtype='poisson'
#
if samplingtype=='poisson':
    print("Poisson-based temporal resampling scheme will be used!")
elif samplingtype=='empirical':
    print("Empirically-based temporal resampling scheme will be used!")
elif samplingtype=='negbinom':
    print("Negative binomial-based temporal resampling scheme will be used! WARNING: this is not tested!")
#
try:
    areatype=cardinfo["POINTAREA"]
    if areatype.lower()=="basin" or areatype.lower()=="watershed":
        areatype='basin'
        try:
            wsmaskshp=str(cardinfo[cardinfo[:,0]=="WATERSHEDSHP",1][0])
            print("You selected 'watershed' for 'POINTAREA', please note that if the watershed is not in a regular lat/lon projection such as EPSG4326/WGS 84, the results will likely be incorrect!")
        except IndexError:
            sys.exit("You specified 'watershed' for 'POINTAREA' but didn't specify 'WATERSHEDSHP'")
    elif areatype.lower()=="point" or areatype.lower()=="grid":
        try:
            ptlat=np.float32(cardinfo[cardinfo[:,0]=="POINTLAT",1][0])
            ptlon=np.float32(cardinfo[cardinfo[:,0]=="POINTLON",1][0])
            ctrlat=ptlat
            ctrlon=ptlon
            areatype="point"
        except ImportError:
            sys.exit("You specified 'point' for 'POINTAREA' but didn't properly specify 'POINTLAT' and 'POINTLON'")
    elif areatype.lower()=="box" or areatype.lower()=="rectangle":
        try:
            box1,box2,box3,box4=cardinfo["BOX_YMIN"],cardinfo["BOX_YMAX"],cardinfo["BOX_XMIN"],cardinfo["BOX_XMAX"]
            boxarea=np.array([box3,box4,box1,box2])
            ctrlat, ctrlon=(box1+box2)/2., (box3+box4)/2
            areatype="box"
        except ImportError:
            sys.exit("You specified 'box' for 'POINTAREA' but didn't properly specify 'BOX_YMIN', 'BOX_YMAX', 'BOX_XMIN', and 'BOX_XMAX'")
    elif areatype.lower()=="pointlist":
        if CreateCatalog:
            sys.exit("POINTLIST is currently not available when creating a new storm catalog.")
        ptlistname=str(cardinfo[cardinfo[:,0]=="POINTLIST",1][0])
        #ptlistname=np.str(sys.argv[3])
        ptlistdat=np.loadtxt(ptlistname,skiprows=1,delimiter=',')
        ptlatlist=ptlistdat[:,0]
        ptlonlist=ptlistdat[:,1]
        npoints_list=len(ptlatlist)
    else:
        sys.exit("unrecognized area type")
except ImportError:
    sys.exit("You didn't specify 'POINTAREA', which needs to either be set to 'watershed', 'point', 'grid', or 'rectangle'")

# transposition
try:
    transpotype=cardinfo["TRANSPOSITION"]
    if transpotype.lower()=='kernel' or transpotype.lower()=='nonuniform':
        transpotype='kernel'
        print("You selected the kernel density-based non-uniform storm transposition scheme!")
    elif transpotype.lower()=='uniform' and areatype.lower()=="pointlist":
        print("You selected the spatially uniform storm transposition scheme and to perform IDF for a list of points!")
    elif transpotype.lower()=='user':
        transpotype='user'
        sys.exit("sadly we aren't set up for the user-supplied transposition scheme yet")
    elif transpotype.lower()=='uniform':
        transpotype='uniform'
        print("You selected the spatially uniform storm transposition scheme!")
    else:
        sys.exit("You specified an unknown resampling method")
except Exception:
    transpotype='uniform'
    print("You didn't specify 'TRANSPOSITION', defaulting to spatially uniform scheme!")

#
# # Use stochastic or deterministic ratio multiplier?
# rescaletype='none'
# try:
#     rescalingtype=cardinfo[cardinfo[:,0]=="ENHANCEDSST",1][0]
#     if rescalingtype.lower()=='stochastic' or rescalingtype.lower()=='deterministic':
#         print("You selected the 'Ratio Rescaling'! This is a more advanced option that should be used with caution")
#         if rescalingtype.lower()=='stochastic':
#             rescaletype='stochastic'
#             print("You selected stochastic ratio rescaling. This has not been thoroughly vetted. Be careful!")
#             if areatype.lower()=='pointlist':
#                 pass
#                 #sys.exit("You selected 'pointlist' for POINTAREA. This is currently not compatible with stochastic rescaling.")
#         if rescalingtype.lower()=='deterministic':
#             rescaletype='deterministic'
#             print("You selected deterministic ratio rescaling. This has not been thoroughly vetted. Be careful!")
#
#         try:
#             rescalingfile=cardinfo[cardinfo[:,0]=="RAINDISTRIBUTIONFILE",1][0]
#             if os.path.isfile(rescalingfile)==False:
#                 sys.exit("The precipitation file specified in 'RAINDISTRIBUTIONFILE' cannot be found!")
#         except IndexError:
#             sys.exit("Even though you 'ratio rescaling', you didn't specify the file of precipitation distributions!")
#     elif rescalingtype.lower()=='dimensionless':
#         print("You selected 'dimensionless SST', modeled after Nathan et al. (2016). This has not been thoroughly vetted. Be careful!")
#         rescaletype='dimensionless'
#         try:
#             rescalingfile=cardinfo[cardinfo[:,0]=="RAINDISTRIBUTIONFILE",1][0]
#             if os.path.isfile(rescalingfile)==False:
#                 sys.exit("The precipitation file specified in 'RAINDISTRIBUTIONFILE' cannot be found!")
#         except IndexError:
#             sys.exit("Even though you 'dimensionless SST', you didn't specify the file of precipitation distributions!")
#     else:
#         print("No rescaling will be performed!")
# except Exception:
#     print("No rescaling will be performed!")
#
#
# try:
#     do_deterministic=cardinfo[cardinfo[:,0]=="MAXTRANSPO",1][0]
#     if do_deterministic.lower():
#         deterministic=True
#     else:
#         deterministic=False
# except Exception:
#     deterministic=False
#
try:
    domain_type=cardinfo["DOMAINTYPE"]
except Exception:
    domain_type='rectangular'

if domain_type.lower()=='rectangular':
    shpdom=False
    ncfdom=False

if domain_type.lower()=='irregular':
    print("Irregular domain type selected!")
    domain_type='irregular'
    ncfdom=False
    shpdom=False
    if CreateCatalog:
        try:
            domainshp=cardinfo["DOMAINSHP"]
            if os.path.isfile(domainshp)==False:
                sys.exit("can't find the transposition domain shapefile!")
            else:
                print("You selected 'irregular' for 'DOMAINTYPE', please note that if the domain shapefile is not in a regular lat/lon projection such as EPSG4326/WGS 84, the results will likely be incorrect!")
                shpdom=True
                ds = shapefile.Reader(domainshp,'rb')
                tempbox= ds.bbox
                inarea=np.array([tempbox[0],tempbox[2],tempbox[1],tempbox[3]],dtype='float32')

        except Exception:
            pass
        try:
            domainncf=cardinfo["DOMAINFILE"]
            if os.path.isfile(domainncf)==False:
                sys.exit("This capability isn't tested. Unclear whether changes are needed due to updating for CF compliance.")
                sys.exit("Can't find the transposition domain NetCDF file!")
            else:
                sys.exit("This capability isn't tested. Unclear whether changes are needed due to updating for CF compliance.")
                print("Domain NetCDF file found!")
                ncfdom=True
                domainmask,domainlat,domainlon=RainyDay.readdomainfile(domainncf)
                res=np.abs(np.mean(domainlat[1:]-domainlat[0:-1]))
                xmin=np.min(np.where(np.sum(domainmask,axis=0)!=0))
                xmax=np.max(np.where(np.sum(domainmask,axis=0)!=0))
                ymin=np.min(np.where(np.sum(domainmask,axis=1)!=0))
                ymax=np.max(np.where(np.sum(domainmask,axis=1)!=0))
                domainmask=domainmask[ymin:(ymax+1),xmin:(xmax+1)]
                inarea=np.array([domainlon[xmin],domainlon[xmax]+res,domainlat[ymax]-res,domainlat[ymin]])
                domainlat=domainlat[ymin:(ymax+1)]
                domainlon=domainlon[xmin:(xmax+1)]
        except Exception:
            pass
    else:
        try:
            domainshp=cardinfo["DOMAINSHP"]
            if os.path.isfile(domainshp)==False:
                sys.exit("can't find the transposition domain shapefile!")
            else:
                print("You selected 'irregular' for 'DOMAINTYPE', please note that if the domain shapefile is not in a regular lat/lon projection such as EPSG4326/WGS 84, the results will likely be incorrect!")
                shpdom=True
        except Exception:
            print("Trouble finding the domain shapefile. Technically we don't need it, so we'll skip this part.")
        yres=np.abs(np.mean(latrange[1:]-latrange[0:-1]))
        xres=np.abs(np.mean(lonrange[1:]-lonrange[0:-1]))
        inarea=np.array([lonrange[0],lonrange[-1]+res,latrange[-1]-res,latrange[0]])

    if ncfdom==False and shpdom==False and CreateCatalog:
        sys.exit("You selected 'irregular' for 'DOMAINTYPE', but didn't provide a shapefile or NetCDF file for the domain!")
else:  ###### Ashar:  Dan I beleive we can remove this condition and add it to the if condition above, this may avoid extra condition.
    print("rectangular domain type selected!")
    domain_type='rectangular'
    if CreateCatalog==False:
        inarea=deepcopy(catarea)
    else:
        try:
            lim1,lim2,lim3,lim4 = cardinfo["Area_Extent"].values()
            inarea=da.array([lim3,lim4,lim1,lim2])
        except ImportError:
            sys.exit("need to specify 'LATITUDE_MIN', 'LATITUDE_MAX', 'LONGITUDE_MIN', 'LONGITUDE_MAX'")


if CreateCatalog==False and da.allclose(inarea,catarea)==False:
    if inarea[0]>=catarea[0] and inarea[1]<=catarea[1] and inarea[2]>=catarea[2] and inarea[3]<=catarea[3]:
        sys.exit("You are changing the domain size, but it fits inside the catalog domain. That's a start, but RainyDay still can't handle it.")
    else:
        sys.exit("You are changing the domain size for an existing catalog. RainyDay can't handle that!")
#
# try:
#     DoDiagnostics=cardinfo[cardinfo[:,0]=="DIAGNOSTICPLOTS",1][0]
#     if DoDiagnostics.lower()=='true':
#         DoDiagnostics=True
#     else:
#         DoDiagnostics=False
# except Exception:
#     DoDiagnostics=True
#
# diagpath=fullpath+'/Diagnostics/'
# if DoDiagnostics==True:
#     if os.path.isdir(wd+scenarioname+'/Diagnostics')==False:
#         #os.system('mkdir %s' %(diagpath))
#         try:
#             os.mkdir(diagpath)
#         except:
#             pass
# #if DoDiagnostics and shpdom:
# #                import geopandas as gpd
# #                dshp = gpd.read_file(domainshp)
# #                #print stateshp1.crs
# #                dshp.crs={}
# #                dshp.to_file(domainshp, driver='ESRI Shapefile')
#
# try:
#     FreqAnalysis=cardinfo[cardinfo[:,0]=="FREQANALYSIS",1][0]
#     if FreqAnalysis.lower()=='true':
#         FreqAnalysis=True
#     else:
#         FreqAnalysis=False
# except Exception:
#     FreqAnalysis=True
#
#
# if areatype.lower()=="pointlist":
#     FreqFile_mean=fullpath+'/'+scenarioname+'_mean.FreqAnalysis'
#     FreqFile_min=fullpath+'/'+scenarioname+'_min.FreqAnalysis'
#     FreqFile_max=fullpath+'/'+scenarioname+'_max.FreqAnalysis'
# else:
#     FreqFile=fullpath+'/'+scenarioname+'_FreqAnalysis.csv'
#
try:
    Scenarios=cardinfo["SCENARIOS"]
    if Scenarios.lower()=='true' and areatype.lower()!="pointlist":
        Scenarios=True
        FreqAnalysis=True
        WriteName=wd+'/'+scenarioname+'/Realizations'
        #os.system('mkdir %s' %(WriteName))
        if os.path.isdir(WriteName)==False:
            os.mkdir(WriteName)
        WriteName=WriteName+'/'+scenarioname
        print("RainyDay will write "+str(nrealizations)+" realizations times "+str(nsimulations)+" years worth of output precipitation scenarios. If this is a big number, this could be very slow!")
    elif Scenarios.lower()=='true' and areatype.lower()=="pointlist":
        print("You specified 'POINTAREA pointlist', but want precipitation scenario outputs. The 'pointlist option' does not support precipitation scenarios!")
        Scenarios=False
    else:
        Scenarios=False
except Exception:
    Scenarios=False
    print("You didn't specify 'SCENARIOS', defaulting to 'false', no scenarios will be written!")
#
#
# # EXLCUDE CERTAIN MONTHS
try:
    exclude=cardinfo["EXCLUDEMONTHS"]
    if type(exclude) == str:
        if exclude.lower() == 'none':

            excludemonths=False
    else:
        excludemonths=exclude
except Exception:
    excludemonths=False
#
# # INCLUDE ONLY CERTAIN YEARS
try:
    includeyr=cardinfo["INCLUDEYEARS"]
    if includeyr.lower()!="all":
        if ',' in includeyr:
            includeyr=includeyr.split(',')
            includeyears=np.empty(len(includeyr),dtype="int32")
            for i in np.arange(0,len(includeyr)):
                includeyears[i]=np.int(includeyr[i])
        elif '-' in includeyr:
            includeyr=includeyr.split('-')
            includeyears=np.arange(np.int(includeyr[0]),np.int(includeyr[1])+1,dtype='int32')
        else:
            includeyears=np.empty((1),dtype="int32")
            includeyears[0]=np.int(includeyr)
    else:
        includeyears=False
except Exception:
    includeyears=False
#
# # MINIMUM RETURN PERIOD THRESHOLD FOR WRITING RAINFALL SCENARIO (THIS WILL HAVE A HUGE IMPACT ON THE SIZE OF THE OUTPUT FILES!)  IN PRINCIPLE IT IS PERHAPS BETTER TO USE AN INTENSITY BUT THIS IS CLEANER!
try:
    RainfallThreshYear=cardinfo["RETURNTHRESHOLD"]
except Exception:
    RainfallThreshYear=1
#
# # DIRECTORY INFO-WHERE DOES THE INPUT DATASET RESIDE?
if CreateCatalog:
    try:
        inpath=cardinfo["RAINPATH"]
        ds = xr.open_mfdataset(inpath, parallel =True)
    except ImportError:
        sys.exit("You didn't specify 'RAINPATH', which is a required field for creating a new storm catalog!")

#
# # SENSITIVITY ANALYSIS OPTIONS:
# try:
#     IntensitySens=cardinfo[cardinfo[:,0]=="SENS_INTENSITY",1][0]
#     if IntensitySens.lower()!='false':
#         IntensitySens=1.+np.float(IntensitySens)/100.
#     else:
#         IntensitySens=1.
# except Exception:
#     IntensitySens=1.
#
# try:
#     FrequencySens=cardinfo[cardinfo[:,0]=="SENS_FREQUENCY",1][0]
#     if FrequencySens.lower()!='false':
#         if samplingtype.lower()!='poisson':
#             sys.exit("Sorry, you can't modify the resampling frequency unless you use the Poisson resampler.")
#         else:
#             FrequencySens=1.0+np.float(FrequencySens)/100.
#     else:
#         FrequencySens=1.
# except Exception:
#     FrequencySens=1.
#     pass
#
# if areatype=='point' or areatype=='pointlist':
#     try:
#         ARFval=cardinfo[cardinfo[:,0]=="POINTADJUST",1][0]
#         if ARFval.lower()!='false':
#             print("Using the POINTADJUST option! Be sure you know what you're doing!")
#             if ',' in ARFval:
#                 ARFval=ARFval.replace(" ","")
#                 arfval=np.array(ARFval.split(','),dtype='float32')
#             else:
#                 arfval=np.array([np.float(ARFval)])
#
#         else:
#             FrequencySens=1.
#     except Exception:
#         arfval=np.array([1.])
#         pass
# else:
#     arfval=np.array([1.])
#
try:
    CalcType=cardinfo["CALCTYPE"]
    if CalcType.lower()=='ams' or CalcType.lower()=='annmax':
        calctype='ams'
        print("You've selected to use Annual Maxima Series.")
    elif CalcType.lower()=='pds' or CalcType.lower()=='partialduration':
        calctype='pds'
        print("You've selected to use Partial Duration Series.")
    else:
        print("Unrecognized entry for CALCTYPE. Defaulting to Annual Maxima Series.")
except Exception:
    calctype='ams'
    print("Nothing provided for CALCTYPE. Defaulting to Annual Maxima Series.")

if calctype=='pds' and Scenarios:
    sys.exit("Unfortunately, RainyDay currently only allows Partial Duration Series for rainfall frequency analysis, not for scenario writing. You can change CALCTYPE to 'AMS' and try again.")
#
#
# # here is where we'll determine if you want to write more than one storm per year to realization files
# # ideally would write the if statement better to provide more guidance if the user provides some weird combination, such as 'CALCTYPE pds' and 'NPERYEAR 5'
if Scenarios and areatype.lower()!="pointlist":
    try:
        nperyear=cardinfo["NPERYEAR"]
        if nperyear!='false' or np.int(nperyear)!=1:
            print("RainyDay will output "+str(np.int(nperyear))+" storms per synthetic year! If this is a big number, this could be very slow!")
            calctype='npyear'
            nperyear=np.int(nperyear)
            if nperyear==1:
                nperyear='false'
    #        else:
    #            print("RainyDay will output "+str(nperyear)+" storms per synthetic year!")
    except Exception:
        pass
#
#
#
try:
    userdistr=cardinfo["INTENSDISTR"]
    if userdistr.lower()=='false':
        userdistr=np.zeros((1),dtype='bool')
    elif len(userdistr.split(','))==3:
        print("reading user-defined precipitation intensity distribution...")
        userdistr=np.array(userdistr.split(','),dtype='float32')
    else:
        sys.exit("There is a problem with the INTENSDISTR entry!")
except Exception:
    userdistr=np.zeros((1),dtype='bool')
    pass

#
#
# if Scenarios:
#     try:
#         if np.any(np.core.defchararray.find(list(cardinfo.keys()),"SPINPERIOD")>-1):
#             pretime=cardinfo[cardinfo[:,0]=="SPINPERIOD",1][0]
#
#             if pretime.lower()=='none' or pretime.lower()=='false':
#                 prependrain=False
#                 print('spinup precipitation will be not included in precipitation scenarios...')
#             else:
#                 try:
#                     pretime=np.int(pretime)
#                     prependrain=True
#                     if pretime==0.:
#                         prependrain=False
#                         print('the spinup time you provided was 0, no spinup precipitation will be included...')
#                     elif pretime<0.:
#                         sys.exit("the spin-up time value is negative.")
#                     elif pretime<1.:
#                         print('the spinup time you provided was short, rounding up to 1 day...')
#                         pretime=1.0
#                     else:
#                         print( np.str(pretime)+' days of spinup precipitation will be included in precipitation scenarios...')
#                 except:
#                     sys.exit('unrecognized value provided for SPINPERIOD.')
#
#         else:
#             prependrain=False
#             print('no SPINPERIOD field was provided.  Spin-up precipitation will be not be included...')
#     except Exception:
#         print('no SPINPERIOD field was provided.  Spin-up precipitation will be not be included...')
#
# try:
#     if np.any(np.core.defchararray.find(list(cardinfo.keys()),"UNCERTAINTY")>-1):
#         spread=cardinfo[cardinfo[:,0]=="UNCERTAINTY",1][0]
#         if spread.lower()=='ensemble':
#             spreadtype='ensemble'
#             print('"ensemble spread" will be calculated for precipitation frequency analysis...')
#         else:
#             try:
#                 int(spread)
#             except:
#                 print('unrecognized value for UNCERTAINTY.')
#             if int(spread)>=0 and int(spread)<=100:
#                 spreadtype='quantile'
#                 quantilecalc=int(spread)
#                 printcalc=(100-quantilecalc)//2
#                 print(str(printcalc)+'th-'+str(quantilecalc+printcalc)+'th interquantile range will be calculated for precipitation frequency analysis...')
#             else:
#                 sys.exit('invalid quantile range for frequency analysis...')
#     else:
#         spreadtype='ensemble'
# except Exception:
#     spreadtype='ensemble'
#
# try:
#     if np.any(np.core.defchararray.find(list(cardinfo.keys()),"RETURNLEVELS")>-1):
#         speclevels=cardinfo[cardinfo[:,0]=="RETURNLEVELS",1][0]
#         if speclevels.lower()=='all':
#             print('using all return levels...')
#             alllevels=True
#         else:
#             alllevels=False
#             if ',' in speclevels:
#                 speclevels=speclevels.split(',')
#                 try:
#                     speclevels=np.float32(speclevels)
#                 except:
#                     print("Non-numeric value provided to RETURNLEVELS.")
#
#                 speclevels=speclevels[speclevels<=nsimulations+0.00001]
#                 speclevels=speclevels[speclevels>=0.99999]
#             else:
#                 sys.exit("The format of RETURNLEVELS isn't recognized.  It should be 'all' or a comma separated list.")
#     else:
#         alllevels=True
# except Exception:
#     alllevels=True
#
#
# try:
#     rotation=False
#     if np.any(np.core.defchararray.find(list(cardinfo.keys()),"ROTATIONANGLE")>-1):
#         try:
#             rotangle=cardinfo[cardinfo[:,0]=="ROTATIONANGLE",1][0]
#         except ImportError:
#             "You're trying to use storm rotation, but didn't specify 'ROTATIONANGLE'"
#         if rotangle.lower()=='none' or rotangle.lower()=='false' or areatype.lower()=="point":
#             rotation=False
#         else:
#             rotation=True
#             if len(rotangle.lower().split(','))!=3:
#                 sys.exit('Unrecognized entry provided for ROTATIONANGLE.  Should be "none" or "-X,+Y,Nangles".')
#             else:
#                 minangle,maxangle,nanglebins=rotangle.split(',')
#                 try:
#                     minangle=np.float32(minangle)
#                 except:
#                     print('Unrecognized entry provided for ROTATIONANGLE minimum angle.')
#                 try:
#                     maxangle=np.float32(maxangle)
#                 except:
#                     print('Unrecognized entry provided for ROTATIONANGLE maximum angle.')
#                 try:
#                     nanglebins=np.int32(nanglebins,dtype="float32")
#                 except:
#                     print('Unrecognized entry provided for ROTATIONANGLE number of bins.')
#
#                 if minangle>0. or maxangle<0.:
#                     sys.exit('The minimum angle should be negative and the maximum angle should be positive.')
# except Exception:
#     rotation=False
#
# if rotation:
#     print("storm rotation will be used...")
#     delarray=[]
#
#
# # should add a "cell-centering" option!
#
# # read in the start and end hour stuff-for doing IDFs for specific times of day-IS THIS ACTUALLY FUNCTIONAL???:
try:
    starthour=cardinfo["STARTHOUR"]
    subday=True
except:
    starthour=0.
#
try:
    endhour=cardinfo["ENDHOUR"]
    subday=True
except:
    endhour=24.
    subday=False

if starthour!=0. and endhour!=24.:
    if (endhour-starthour)<duration:
        duration=endhour-starthour
        sys.exit("Not sure I ever finished writing this capability!")

try:
    if np.any(np.core.defchararray.find(list(cardinfo.keys()),"DURATIONCORRECTION")>-1):
        durcorr=cardinfo["DURATIONCORRECTION"]
        if durcorr.lower()=='false':
            durcorrection=False
            print('DURATIONCORRECTION set to "false". No correction will be used.')
        elif durcorr.lower()=='true':
            durcorrection=True
        else:
            print('Invalid option provided for DURATIONCORRECTION (should be "true" or "false"). Defaulting to "false"!')
            durcorrection=False
    else:
        print('No value was given for DURATIONCORRECTION. If applicable, it will default to "false"!')
        durcorrection=False
except Exception:
    durcorrection=False

# Added by DBW, 16 Feb 2021, specifically to support Dr. Emad Habib's team in performing ARF analyses.
# If ARFANALYSIS is used, DURATIONCORRECTION will be turned on automatically, and the scenarios will have a duration equal to DURATION.
#This is sensible if the objective is to create rainfall scenario files for subsequent analysis (ARFs, scaling properties, etc.).
# If ARFANALYSIS is turned off and DURATIONCORRECTION is turned on, then the rainfall scenarios will under most circumstances have a longer duration than DURATION.
# This latter approach (DURATIONCORRECTION on, ARFANALYSIS off) is a more justifiable way of generating rainfall scenarios for flood modeling.
# Also, this ARFANALYSIS option isn't set up to work with 'pointlists'.
try:
    if np.any(np.core.defchararray.find(list(cardinfo.keys()),"ARFANALYSIS")>-1):
        arfcorr=cardinfo[cardinfo[:,0]=="ARFANALYSIS",1][0]
        if arfcorr.lower()=='false':
            arfcorrection=False
        elif arfcorr.lower()=='true' and areatype.lower()!="pointlist":
            durcorrection=True
            arfcorrection=True
            print("ARFANALYSIS set to 'true'. DURATIONCORRECTION will be used. This is only really recommended in you want to create rainfall scenarios for ARF or other spatial analysis. Make sure you know what you're doing!")
        elif prependrain:
            sys.exit("You've selected to do ARFANALYSIS and to pre-pend rainfall. This suggests that you don't know what you're doing.")
        else:
            print('Invalid option provided for ARFANALYSIS (should be "true" or "false"). Defaulting to "false"!')
            arfcorrection=False
    else:
        arfcorrection=False
except Exception:
    arfcorrection=False

if CreateCatalog and durcorrection:
    catduration=max([72.,3.*duration])    # this will be slow!
    print('Since DURCORRECTION will be used and a storm catalog will be created, the duration of the catalog will be '+"{0:0.2f}".format(catduration)+' hours')
elif CreateCatalog:
    catduration=duration


#sys.exit('decide what to do with duration correction and timeseparation!')
try:
    timeseparation=cardinfo["TIMESEPARATION"]
except Exception:
    timeseparation=0.

islassiter=False
try:
    if np.any(np.core.defchararray.find(list(cardinfo.keys()),"ISLASSITER")>-1):
        islassiter=cardinfo[cardinfo[:,0]=="ISLASSITER",1][0]
        if islassiter.lower()!="false":
            islassiter=True
            print("Using Lassiter-style Netcdf files!")
        else:
            islassiter=False
except Exception:
    pass
try:
    if np.any(np.core.defchararray.find(list(cardinfo.keys()),"ISFITZGERALD")>-1):
        islassiter=cardinfo[cardinfo[:,0]=="ISFITZGERALD",1][0]
        if islassiter.lower()!="false":
            islassiter=True
            print("Using FitzGerald-style Netcdf files!")
        else:
            islassiter=False
except Exception:
    pass


#==============================================================================
# THIS BLOCK CONFIGURES SEVERAL THINGS
#==============================================================================

initseed=0
da.random.seed(initseed)
global rainprop
rainprop=deepcopy(emptyprop)
#%%

#==============================================================================
# CREATE NEW STORM CATALOG, IF DESIRED
#==============================================================================
if CreateCatalog:
    print("creating a new storm catalog...")

    ds,nyears=RainyDay.createfilelist(ds, includeyears, excludemonths)
    if defaultstorms:
        nstorms=nyears*20


    # GET SUBDIMENSIONS, ETC. FROM THE NETCDF FILE RATHER THAN FROM RAINPROPERTIES
    rainprop.spatialres,rainprop.dimensions,rainprop.bndbox,rainprop.timeres,rainprop.nodata=RainyDay.rainprop_setup(ds,lassiterfile=islassiter)
    spatres=rainprop.spatialres[0].values


    #==============================================================================
    # SET UP THE SUBGRID INFO
    #==============================================================================
    # 'subgrid' defines the transposition domain

    rainprop.subextent,rainprop.subind,rainprop.subdimensions=RainyDay.findsubbox(inarea,rainprop)
    if ncfdom and (da.any(domainmask.shape!=rainprop.subdimensions)):
        sys.exit("Something went terribly wrong :(")        # this shouldn't happen



#==============================================================================
# IF A STORM CATALOG ALREADY EXISTS, USE IT
#==============================================================================

else:
    print("Using the existing storm catalog...")

    rainprop.spatialres=[xres,yres]
    if len(cattime[-1,1:])!=0:
        rainprop.timeres=np.float32(np.mean(cattime[-1,1:]-cattime[-1,0:-1]))
    else:
        rainprop.timeres=timeres

    rainprop.nodata=np.unique(catrain[catrain<0.])

    delt=np.timedelta64(cattime[0,-1]-(cattime[0,0]))
    catduration=(delt.astype('int')+rainprop.timeres)/60.
    if defaultstorms==False:
        if nstorms>cattime.shape[0]:
            print("WARNING: The storm catalog has fewer storms than the specified nstorms")
            nstorms=cattime.shape[0]
    else:
        darr=np.array(pd.to_datetime(cattime.ravel()).year)
        nstorms_default=(np.max(darr)-np.min(darr)+1)*25
        if nstorms_default>cattime.shape[0]:
            print("WARNING: The storm catalog has fewer storms than the default number of storms")
            nstorms=cattime.shape[0]
        else:
            nstorms=nstorms_default


    if ncfdom and (np.any(domainmask.shape!=catrain.shape[2:])):
        sys.exit("The domain mask and the storm catalog are different sizes. Do you know what you're doing?")

    rainprop.bndbox=catarea
    rainprop.subextent=rainprop.bndboxreadrealization

    #rainprop.subextent,rainprop.subind,rainprop.subdimensions=RainyDay.findsubbox(inarea,rainprop)

    rainprop.subind=np.array([0,len(lonrange)-1,len(latrange)-1,0],dtype='int32')
    rainprop.dimensions=np.array([len(latrange),len(lonrange)],dtype='int32')

    #catrain=catrain[:,:,rainprop.subind[3]:rainprop.subind[2]+1,rainprop.subind[0]:rainprop.subind[1]+1]
    #latrange=latrange[rainprop.subind[0]:rainprop.subind[1]+1]
    #lonrange=lonrange[rainprop.subind[3]:rainprop.subind[2]+1]

    rainprop.subdimensions=rainprop.dimensions

if int(duration*60/rainprop.timeres)<=0:
    sys.exit("it appears that you specified a duration shorter than the temporal resolution of the input data!")


if timeseparation<=0. and durcorrection==False:
    timeseparation=duration
elif timeseparation>0. and durcorrection==False:
    timeseparation=timeseparation+duration
elif durcorrection:
    timeseparation=max([timeseparation+duration,catduration])     ### This is

# spatres=rainprop.spatialres[0]
ingridx,ingridy=np.meshgrid(da.arange(rainprop.subextent[0],rainprop.subextent[1]-spatres/1000.,spatres),np.arange(rainprop.subextent[3],rainprop.subextent[2]+spatres/1000.,-spatres))
lonrange=ingridx[0,:]
latrange=ingridy[:,0]

#============================================================================
# Do the setup to run for specific times of day!    ### Won't work with this form of xarray (have to change the code a little bit
#=============================================================================

if CreateCatalog==False:
    tdates = pd.DatetimeIndex(cattime[:,0])
    tdates.year
    nyears=np.max(np.array(tdates.year))-np.min(np.array(tdates.year))+1
if starthour==0 and endhour==24:
    hourinclude=da.ones((int(24*60/rainprop.timeres)),dtype='int32')
else:
    sys.exit("Restrictions to certain hours isn't currently tested or supported")
    try:
        sys.exit("need to fix this")
        _,temptime,_,_=RainyDay.readnetcdf(flist[0],inbounds=rainprop.subind,lassiterfile=islassiter)
    except Exception:
        sys.exit("Can't find the input files necessary for setup to calculate time-of-day-specific IDF curves")

    starthour=starthour+np.float(rainprop.timeres/60)   # because the netcdf file timestamps correspond to the end of the accumulation period
    hourinclude=np.zeros((24*60/rainprop.timeres),dtype='int32')
    temphour=np.zeros((24*60/rainprop.timeres),dtype='float32')

    for i in np.arange(0,len(temptime)):
        temphour[i]=temptime[i].astype(object).hour

    # the following line will need to be adapted, due to UTC vs. local issues
    hourinclude[np.logical_and(np.greater_equal(temphour,starthour),np.less_equal(temphour,endhour))]=1
    if len(hourinclude)!=len(temptime):
        sys.exit("Something is wrong in the hour exclusion calculation!")

hourinclude=hourinclude.astype('bool')

#temptime[hourinclude] # checked, seems to be working right


#==============================================================================
# SET UP GRID MASK
#==============================================================================

print("setting up the grid information and masks...")
if areatype.lower()=="basin" or areatype.lower()=="watershed":
    if os.path.isfile(wsmaskshp)==False:
        sys.exit("can't find the basin shapefile!")
    else:
        catmask=RainyDay.rastermask(wsmaskshp,rainprop,'fraction')
        #catmask=catmask.reshape(ingridx.shape,order='F')

elif areatype.lower()=="point":
    catmask=np.zeros((rainprop.subdimensions))
    yind=np.where((latrange-ptlat)>0)[0][-1]
    xind=np.where((ptlon-lonrange)>0)[0][-1]
    if xind==0 or yind==0:
        sys.exit('the point you defined is too close to the edge of the box you defined!')
    else:
        catmask[yind,xind]=1.0

elif areatype.lower()=="pointlist":
    catmask=np.zeros((rainprop.subdimensions))
    yind_list=np.zeros_like(ptlatlist,dtype='int32')
    xind_list=np.zeros_like(ptlatlist,dtype='int32')
    for pt in np.arange(0,npoints_list):
        yind_list[pt]=np.where((latrange-ptlatlist[pt])>0)[0][-1]
        xind_list[pt]=np.where((ptlonlist[pt]-lonrange)>0)[0][-1]
    if np.any(yind_list==0) or np.any(xind_list==0):
        sys.exit('the point you defined is too close to the edge of the box you defined!')
    else:
        catmask[yind_list[1],xind_list[1]]=1.0          # the idea of the catmask gets a little goofy in this situation

elif areatype.lower()=="box":
    print("I am here at areatype")
    x,y = rainprop.spatialres[0].values, rainprop.spatialres[1].values
    finelat = np.arange(latrange[0],latrange[-1]- y+ x/1000,- y/25)
    finelon = np.arange(lonrange[0],lonrange[-1]+ x - x/1000, x/25)
    subindy= np.logical_and(finelat>boxarea[2]+ y/1000,finelat<boxarea[3]+ y/1000)
    subindx= np.logical_and(finelon>boxarea[0]-x/1000,finelon<boxarea[1]-x/1000)

    tx,ty=np.meshgrid(subindx,subindy)
    catmask=np.array(da.logical_and(tx==True,ty==True),dtype='float32')

    if len(finelat[subindy])<25 and len(finelon[subindx])<25:
        print('WARNING: you set POINTAREA to "box", but the box is smaller than a single pixel.  This is not advised.  Either set POINTAREA to "point" or increase the size of the box.'   )

    if len(finelat[subindy])==1 and len(finelon[subindx])==1:
        catmask=np.zeros((rainprop.subdimensions))
        yind=np.where((latrange-ptlat)>0)[0][-1]
        xind=np.where((ptlon-lonrange)>0)[0][-1]
        if xind==0 or yind==0:
            sys.exit('the point you defined is too close to the edge of the box you defined!')
        else:
            catmask[yind,xind]=1.0
    else:
        def block_mean(ar, fact):
            assert isinstance(fact, int), type(fact)
            sx, sy = ar.shape
            X, Y = np.ogrid[0:sx, 0:sy]
            regions = sy//fact * (X//fact) + Y//fact
            res = ndimage.mean(ar, labels=regions, index=np.arange(regions.max() + 1))
            res.shape = (sx//fact, sy//fact)
            return res

        catmask=block_mean(catmask,25)
    print("I am at the end of areatype")
    # this scheme is a bit of a numerical approximation but I doubt it makes much practical difference
else:
    sys.exit("unrecognized area type!")


# TRIM THE GRID DOWN
print("I am trimming")
csum, rsum=da.where(da.sum(catmask,axis=0)==0), da.where(da.sum(catmask,axis=1)==0)

xmin=da.where(da.sum(catmask,axis=0)!=0)[0].min()
xmax=da.where(da.sum(catmask,axis=0)!=0)[0].max()
ymin=da.where(da.sum(catmask,axis=1)!=0)[0].min()
ymax=da.where(da.sum(catmask,axis=1)!=0)[0].max()

trimmask  =da.delete(catmask,csum,axis=1)
trimmask  =da.delete(trimmask,rsum,axis=0)
maskwidth =trimmask.shape[1]
maskheight=trimmask.shape[0]
trimmask  =da.array(trimmask,dtype='float32')
catmask   =da.array(catmask,dtype='float32')

timeseparation=np.timedelta64(np.int(timeseparation*60.),'m')
timestep      =np.timedelta64(np.int(rainprop.timeres),'m')

mnorm=da.sum(trimmask)

xlen=rainprop.subdimensions[1]-maskwidth+1
ylen=rainprop.subdimensions[0]-maskheight+1
print("Trimming ends")

if ncfdom:
    ingridx_dom,ingridy_dom=np.meshgrid(domainlon,domainlat)

    ingrid_domain=np.column_stack((ingridx_dom.flatten(),ingridy_dom.flatten()))
    #ingridx,ingridy=np.meshgrid(np.arange(rainprop.subextent[0],rainprop.subextent[1]-rainprop.spatialres[0]/1000,rainprop.spatialres[0]),np.arange(rainprop.subextent[3],rainprop.subextent[2]+rainprop.spatialres[1]/1000,-rainprop.spatialres[1]))
    grid_rainfall=np.column_stack((ingridx.flatten(),ingridy.flatten()))

    delaunay=sp.spatial.qhull.Delaunay(ingrid_domain)
    interp=sp.interpolate.NearestNDInterpolator(delaunay,domainmask.flatten())

    # in case you want it back in a rectangular grid:
    domainmask=da.reshape(interp(grid_rainfall),ingridx.shape)

elif domain_type.lower()=='irregular' and shpdom and CreateCatalog:
    domainmask=RainyDay.rastermask(domainshp,rainprop,'simple').astype('float32')
    #domainmask=catmask.reshape(ingridx.shape,order='F')
elif domain_type.lower()=='rectangular':
    domainmask=da.ones((catmask.shape),dtype='float32')
else:
    pass


if catmask.shape!=domainmask.shape:
    sys.exit("Oh dear, 'catmask' and 'domainmask' aren't the same size!")

if np.allclose(np.multiply(catmask,domainmask),0.) and areatype.lower()!="pointlist":
    sys.exit("it looks as if the location specified in 'POINTAREA' is outside of the transposition domain!")

# exclude points that are outside of the transposition domain:
if areatype.lower()=="pointlist" and domain_type.lower()=='irregular':
    keeppoints=np.ones_like(yind_list,dtype='bool')
    for pt in np.arange(0,npoints_list):
        if domainmask[yind_list[pt],xind_list[pt]]==0:
            keeppoints[pt]=False
elif areatype.lower()=="pointlist" and domain_type.lower()=='rectangular':
    keeppoints=np.ones_like(yind_list,dtype='bool')
    for pt in np.arange(0,npoints_list):
        if ptlatlist[pt]>rainprop.subextent[3] or ptlatlist[pt]<rainprop.subextent[2] or ptlonlist[pt]>rainprop.subextent[1] or ptlonlist[pt]<rainprop.subextent[0]:
            keeppoints[pt]=False


if areatype.lower()=="pointlist":
    if np.any(keeppoints)==False:
        sys.exit("it looks as if none of the points in specified in 'POINTLIST' are inside the transposition domain!")
    else:
        yind_list=yind_list[keeppoints]
        xind_list=xind_list[keeppoints]
        ptlatlist=ptlatlist[keeppoints]
        ptlonlist=ptlonlist[keeppoints]
        npoints_list=len(yind_list)



#################################################################################
# STEP 1: CREATE STORM CATALOG
#################################################################################

if CreateCatalog:
    print("reading precipitation files...")

    #==============================================================================
    # SET UP OUTPUT VARIABLE
    #==============================================================================

    rainarray=np.zeros((int(catduration*60/rainprop.timeres),rainprop.subdimensions[0],rainprop.subdimensions[1]),dtype='float32')

    rainsum=da.zeros((ylen,xlen),dtype='float32'); rainarray[:]=np.nan

    raintime=np.empty((int(catduration*60/rainprop.timeres)),dtype='datetime64[m]')
    raintime[:]=np.datetime64(datetime(1900,1,1,0,0,0))

    catmax=da.zeros((nstorms),dtype='float32')

    catrain=da.zeros((nstorms,int(catduration*60/rainprop.timeres),rainprop.subdimensions[0],rainprop.subdimensions[1]),dtype='float32')
    cattime=np.empty((nstorms,int(catduration*60/rainprop.timeres)),dtype='datetime64[m]')
    cattime[:]=np.datetime64(datetime(1900,1,1,0,0,0))
    catloc=da.empty((nstorms),dtype='float32')

    catx,caty=da.zeros((nstorms),dtype='int32'),da.zeros((nstorms),dtype='int32')

    #==============================================================================
    # READ IN RAINFALL
    #==============================================================================

    start = time.time()
    inrain,intime,_,_= RainyDay.readnetcdf(ds,inbounds=rainprop.subind,lassiterfile=islassiter)
    inrain=inrain.where(inrain>0)
    # intime=intime[hourinclude]
    # inrain[inrain<0.]=np.nan
    # print('Processing file '+str(i+1)+' out of '+str(len(flist))+' ('+"{0:0.0f}".format(100*(i+1)/len(flist))+'%): '+infile.split('/')[-1])
    # THIS FIRST PART BUILDS THE STORM CATALOG
    print(" Analysing and creating storms")
    for k in range(0,len(intime)):
        starttime      =intime[k]-np.timedelta64(int(catduration*60.),'m')
        raintime[-1]   =intime[k]
        rainarray[-1,:]=inrain[k,:]
        #rainarray[-1,:]=np.reshape(inrain[k,:],(rainprop.subdimensions[0],rainprop.subdimensions[1]))
        subtimeind=da.where(np.logical_and(raintime>starttime,raintime<=raintime[-1]))
        subtime   =np.arange(raintime[-1],starttime,-timestep)[::-1]
        temparray =da.squeeze(da.nansum(rainarray[subtimeind,:],axis=1))

        if da.any(da.greater(temparray,np.min(catmax))): # DBW-added this if statement on 10112022. It seems like this should speed things up!
            if domain_type=='irregular':
                rainmax,ycat,xcat=RainyDay.catalogNumba_irregular(temparray,trimmask,xlen,ylen,maskheight,maskwidth,rainsum,domainmask)
            else:
                rainmax,ycat,xcat=RainyDay.catalogNumba(temparray,trimmask,xlen,ylen,maskheight,maskwidth,rainsum)

            minind=np.argmin(catmax).compute()
            tempmin=catmax[minind]
            if rainmax>tempmin:
                checksep=intime[k]-cattime[:,-1]
                if (checksep<timeseparation).any():
                    checkind = da.where(checksep<timeseparation)[0].compute()[0]
                    if rainmax>=catmax[checkind]:
                        catmax[checkind]   =rainmax
                        cattime[checkind,:]=subtime
                        catx[checkind]     =xcat
                        caty[checkind]     =ycat
                        catrain[checkind,:]=rainarray
                else:
                    catmax[minind]   =rainmax
                    cattime[minind,:]=subtime
                    catx[minind]     =xcat
                    caty[minind]     =ycat
                    catrain[minind,:]=rainarray


        rainarray[0:-1,:]=rainarray[1:int(catduration*60/rainprop.timeres),:]
        raintime[0:-1]=raintime[1:int(catduration*60/rainprop.timeres)]

    sind=np.argsort(catmax)
    cattime=cattime[sind,:]
    catx=catx[sind]
    caty=caty[sind]
    catrain=catrain[sind,:]
    catmax=catmax[sind]/mnorm*rainprop.timeres/60.

    end = time.time()
    print("catalog timer: "+"{0:0.2f}".format((end - start)/60.)+" minutes")

    # WRITE CATALOG
    print("writing storm catalog...")
    RainyDay.writecatalog(scenarioname,catrain,catmax,catx,caty,cattime,latrange,lonrange,catalogname,nstorms,catmask,parameterfile,domainmask,timeresolution=rainprop.timeres)

