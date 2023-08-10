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
#%%
if FreqAnalysis:
    print("resampling and transposing...")
    
    if np.all(includeyears==False):
        nyears=len(range(min(cattime[:,-1].astype('datetime64[Y]').astype(int)),max(cattime[:,-1].astype('datetime64[Y]').astype(int))+1))
    else:
        nyears=len(includeyears)
        
        
    # resampling counts options:
    if samplingtype.lower()=='poisson':
        lrate=len(catmax)/nyears*FrequencySens                  
        ncounts=np.random.poisson(lrate,(nsimulations,nrealizations))
        cntr=0
        #ncounts[ncounts==0]=1
        if calctype.lower()=='npyear' and lrate<nperyear:   
            sys.exit("You specified to write multiple storms per year, but you specified a number that is too large relative to the resampling rate!")
    elif samplingtype.lower()=='negbinom':
        sys.exit("Sorry, the negative binomial resampling isn't set up yet :(")
        _,yrscount=np.unique(cattime[:,-1].astype('datetime64[Y]').astype(int)+1970,return_counts=True)
        if len(yrscount)<nyears:
            yrscount=np.append(yrscount,np.ones(nyears-len(yrscount),dtype='int32'))
        rparam=np.mean(yrscount)*np.mean(yrscount)/(np.var(yrscount)-np.mean(yrscount))   
        pparam=1.-np.mean(yrscount)/np.var(yrscount)
        ncounts=np.random.negative_binomial(rparam,pparam,size=(nsimulations,nrealizations)) 
        if calctype.lower()=='npyear' and np.mean(yrscount)<nperyear :   
            sys.exit("You specified to write multiple storms per year, but you specified a number that is too large relative to the resampling rate!")
             	
    else:
        _,yrscount=np.unique(cattime[:,-1].astype('datetime64[Y]').astype(int)+1970,return_counts=True)
        if len(yrscount)<nyears:
            yrscount=np.append(yrscount,np.ones(nyears-len(yrscount),dtype='int32'))
        ncounts=np.random.choice(yrscount,(nsimulations,nrealizations),replace=True)   
        if calctype.lower()=='npyear' and np.mean(yrscount)<nperyear :   
            sys.exit("You specified to write multiple storms per year, but you specified a number that is too large relative to the resampling rate!")
        #ncounts[ncounts==0]=1
            
    whichstorms=np.empty((np.nanmax(ncounts),ncounts.shape[0],ncounts.shape[1]),dtype='int32')
    whichstorms[:]=-9999
    
    if rotation==True:
        randangle=(maxangle-minangle)*np.random.random_sample(((np.nanmax(ncounts),ncounts.shape[0],ncounts.shape[1])))+minangle
        
        angbins=np.linspace(minangle,maxangle,nanglebins)
        angs=math.pi/180*angbins
        anglebins=np.digitize(randangle.ravel(),angbins).reshape(np.nanmax(ncounts),ncounts.shape[0],ncounts.shape[1])
    
    
    # DOES THIS PROPERLY HANDLE STORM EXCLUSIONS???  I think so...
    for i in range(0,np.nanmax(ncounts)):
        whichstorms[i,ncounts>=i+1]=np.random.randint(0,nstorms,(len(ncounts[ncounts>=i+1]))) # why was this previously "nstorms-1"??? Bug?
    
    
    # the next three lines were commented out when adding the "pointlist" option
    #whichrain=np.zeros((whichstorms.shape),dtype='float32')
    #whichx=np.zeros((whichstorms.shape),dtype='int32')
    #whichy=np.zeros((whichstorms.shape),dtype='int32')
    
    if areatype.lower()=="pointlist":
        whichx=np.zeros((whichstorms.shape[0],whichstorms.shape[1],whichstorms.shape[2],npoints_list),dtype='int32')
        whichy=np.zeros((whichstorms.shape[0],whichstorms.shape[1],whichstorms.shape[2],npoints_list),dtype='int32')
        whichrain=np.zeros((whichstorms.shape[0],whichstorms.shape[1],whichstorms.shape[2],npoints_list),dtype='float32')
    else:
        whichx=np.zeros((whichstorms.shape[0],whichstorms.shape[1],whichstorms.shape[2],1),dtype='int32')
        whichy=np.zeros((whichstorms.shape[0],whichstorms.shape[1],whichstorms.shape[2],1),dtype='int32')  
        whichrain=np.zeros((whichstorms.shape[0],whichstorms.shape[1],whichstorms.shape[2],1),dtype='float32')
        whichstep=np.zeros((whichstorms.shape[0],whichstorms.shape[1],whichstorms.shape[2],1),dtype='int32')
        
    if durcorrection:
        whichtimeind=np.zeros((whichstorms.shape),dtype='float32')
    
    
    if transpotype=='uniform' and domain_type=='irregular':
        domainmask[-maskheight:,:]=0.
        domainmask[:,-maskwidth:]=0.
        xmask,ymask=np.meshgrid(np.arange(0,domainmask.shape[1],1),np.arange(0,domainmask.shape[0],1))
        xmask=xmask[np.equal(domainmask,True)]
        ymask=ymask[np.equal(domainmask,True)]
        
    if rescaletype=='stochastic' or rescaletype=='deterministic' or rescaletype=='dimensionless':
        whichmultiplier=np.empty_like(whichrain)
        whichmultiplier[:]=np.nan

        
    #==============================================================================
    # If you're using intensity-dependent resampling, get ready for it!
    #==============================================================================


    if rescaletype=='stochastic' or rescaletype=='deterministic':
        smoothsig=5
        
        print("reading in precipitation intensity data...")
        intenserain,_,intenselat,intenselon=RainyDay.readintensityfile(rescalingfile)
        intensemask=np.equal(np.sum(intenserain,axis=0),0.)
        intenserain[:,intensemask]=np.nan
        int_xmin=np.abs(intenselon-rainprop.bndbox[0]).argmin()
        int_ymin=np.abs(intenselat-rainprop.bndbox[3]).argmin()
        int_xmax=np.abs(intenselon-rainprop.bndbox[1]).argmin()
        int_ymax=np.abs(intenselat-rainprop.bndbox[2]).argmin()
        #int_xmax=np.abs(intenselon-rainprop.bndbox[1]).argmin()+1
        #int_ymax=np.abs(intenselat-rainprop.bndbox[2]).argmin()+1
        intenserain=intenserain[:,int_ymin:int_ymax,int_xmin:int_xmax]
        intenselat=intenselat[int_ymin:int_ymax]
        intenselon=intenselon[int_xmin:int_xmax]
        
        #intenserain[np.equal(intenserain,0.)]=np.nan
        nintstorms=np.min((intenserain.shape[0],2*nyears))
        intenserain=intenserain[-nintstorms:,:]

        if np.array_equal(intenselat,latrange)==False or np.array_equal(intenselon,lonrange)==False:  
            intensegridx,intensengridy=np.meshgrid(intenselon,intenselat)        
            ingrid_intense=np.column_stack((intensegridx.flatten(),intensengridy.flatten())) 
            grid_out=np.column_stack((ingridx.flatten(),ingridy.flatten())) 
            delaunay=sp.spatial.qhull.Delaunay(ingrid_intense)
            tempintense=np.empty((intenserain.shape[0],ingridx.shape[0],ingridx.shape[1]),dtype='float32')
            for i in range(0,intenserain.shape[0]):
                interp=sp.interpolate.LinearNDInterpolator(delaunay,intenserain[i,:].flatten(),fill_value=np.nan)
                tempintense[i,:]=np.reshape(interp(grid_out),ingridx.shape)
            intenserain=tempintense 
        
        #intenserain[:,np.equal(domainmask,0.)]=np.nan
      

        # the stochastic multiplier approach uses the log of the rainfall:
        intenserain=np.log(intenserain)
        hometemp=np.nansum(np.multiply(intenserain,catmask),axis=(1,2))/mnorm
        xlen_wmask=rainprop.subdimensions[1]-maskwidth+1
        ylen_wmask=rainprop.subdimensions[0]-maskheight+1
                
        
        # there is some goofy numba stuff below...
        if maskheight>1 or maskwidth>1:
            tempintense=np.empty((intenserain.shape[0],ylen_wmask,xlen_wmask),dtype='float32')
            intenserain=RainyDay.intenseloop(intenserain,tempintense,xlen_wmask,ylen_wmask,maskheight,maskwidth,trimmask,mnorm,domainmask)
        intensecorr=np.empty((ylen_wmask,xlen_wmask),dtype='float32')    
        intensecorr=RainyDay.intense_corrloop(intenserain,intensecorr,hometemp,xlen_wmask,ylen_wmask,mnorm,domainmask)
            
        #homemean=np.mean(homerain)
        #homestd=np.std(homerain)
        intensemean=np.mean(intenserain,axis=0)
        intensestd=np.std(intenserain,axis=0) 
        
        intensemean=RainyDay.mysmoother(intensemean,sigma=[smoothsig,smoothsig])
        intensestd=RainyDay.mysmoother(intensestd,sigma=[smoothsig,smoothsig])
        intensecorr=RainyDay.mysmoother(intensecorr,sigma=[smoothsig,smoothsig])
        
        if maskheight>1 or maskwidth>1:
            homemean=np.nansum(np.multiply(intensemean,catmask[:-maskheight+1,:-maskwidth+1]),axis=(0,1))/mnorm
            homestd=np.nansum(np.multiply(intensestd,catmask[:-maskheight+1,:-maskwidth+1]),axis=(0,1))/mnorm        
        else:    
            homemean=np.nansum(np.multiply(intensemean,catmask),axis=(0,1))/mnorm
            homestd=np.nansum(np.multiply(intensestd,catmask),axis=(0,1))/mnorm
        
        # just in case you don't have any data to inform the rescaling:
        intensemean[np.isneginf(intensemean)]=homemean
        intensestd[np.isneginf(intensestd)]=0.
        intensecorr[np.isneginf(intensecorr)]=1.0
        
        intensemean[np.isnan(intensemean)]=homemean
        intensestd[np.isnan(intensestd)]=0.
        intensecorr[np.isnan(intensecorr)]=1.0
        
    elif rescaletype=='dimensionless':
        print("reading in precipitation map for dimensionless SST...")
        if '.nc' in rescalingfile:
            sys.exit('need to set this up')
            #intenserain,_,intenselat,intenselon=RainyDay.readintensityfile(rescalingfile)
        elif '.asc' in rescalingfile:
            asciigrid,ncols,nrows,xllcorner,yllcorner,cellsize=RainyDay.read_arcascii(rescalingfile)
            dlsstarea=[xllcorner,xllcorner+ncols*cellsize,yllcorner,yllcorner+nrows*cellsize]
            atlasgridx,atlasgridy=np.meshgrid(np.arange(dlsstarea[0],dlsstarea[1]-cellsize/10.,cellsize),np.arange(dlsstarea[3],dlsstarea[2]+cellsize/10.,-cellsize))
            atlas14_domain=np.column_stack((atlasgridx.flatten(),atlasgridy.flatten())) 
        
            delaunay=sp.spatial.qhull.Delaunay(atlas14_domain)
            interp=sp.interpolate.LinearNDInterpolator(delaunay,asciigrid.flatten(),fill_value=np.nan)
            
            grid_out=np.column_stack((ingridx.flatten(),ingridy.flatten())) 
            atlas_regridded=np.reshape(interp(grid_out),ingridx.shape) 
            atlas_regridded=np.log(atlas_regridded)
            if areatype.lower()!='pointlist':
                hometemp=np.nansum(np.multiply(atlas_regridded,catmask))/mnorm
            else:
                hometemp=np.nanmean(atlas_regridded[domainmask==True])
            atlas_regridded[np.isnan(atlas_regridded)]=hometemp
   
        else:
            sys.exit('Unrecognized file format provided for dimensionless SST')   
    else:
        pass
        
   
    # here is the main resampling and transposition loop
    for i in np.arange(0,nstorms):
        print('Resampling and transposing storm '+str(i+1)+' out of '+str(nstorms)+' ('"{0:0.0f}".format(100*(i+1)/nstorms)+'%)')
        # UNIFORM RESAMPLING
        if transpotype=='uniform' and domain_type=='rectangular':
            whichx[whichstorms==i,0]=np.random.randint(0,np.int(rainprop.subdimensions[1])-maskwidth+1,len(whichx[whichstorms==i]))
            whichy[whichstorms==i,0]=np.random.randint(0,np.int(rainprop.subdimensions[0])-maskheight+1,len(whichy[whichstorms==i]))
     
        # KERNEL-BASED AND INTENSITY-BASED RESAMPLING (ALSO NEEDED FOR IRREGULAR TRANSPOSITION DOMAINS)
        elif transpotype=='kernel':
            rndloc=np.array(np.random.random_sample(len(whichx[whichstorms==i])),dtype='float32')
            tempx=np.empty((len(rndloc)),dtype='int32')
            tempy=np.empty((len(rndloc)),dtype='int32')
            for pt in np.arange(0,whichx.shape[3]):
                whichx[whichstorms==i,pt],whichy[whichstorms==i,pt]=RainyDay.numbakernel_fast(rndloc,cumkernel[:,:,pt],tempx,tempy,rainprop.subdimensions[1])

        if transpotype=='uniform' and domain_type=='irregular':
            rndloc=np.random.randint(0,np.sum(np.equal(domainmask,True)),np.sum(whichstorms==i))
            for pt in np.arange(0,whichx.shape[3]):
                whichx[whichstorms==i,pt]=xmask[rndloc].reshape(len(xmask[rndloc]))
                whichy[whichstorms==i,pt]=ymask[rndloc].reshape(len(ymask[rndloc]))
        
        # SET UP MANUAL PDF RESAMPLING
        elif transpotype=='manual':  
            sys.exit("not configured for manually supplied pdf yet!")
    
        if durcorrection:
            passrain=np.array(RainyDay.rolling_sum(catrain[i,:], int(duration*60/rainprop.timeres)),dtype='float32')
            
        else:
            passrain=np.nansum(catrain[i,:],axis=0)         # time-average the rainfall
    
        if rotation: 
            print('rotating storms for transposition, '+str(100*(i+1)/nstorms)+'% complete...')
            delarray.append([])
             
            xctr=catx[i]+maskwidth/2.
            yctr=caty[i]+maskheight/2.
            xlinsp=np.linspace(-xctr,rainprop.subdimensions[1]-xctr,rainprop.subdimensions[1])
            ylinsp=np.linspace(-yctr,rainprop.subdimensions[0]-yctr,rainprop.subdimensions[0])
            ingridx,ingridy=np.meshgrid(xlinsp,ylinsp)
            ingridx=ingridx.flatten()
            ingridy=ingridy.flatten()
            outgrid=np.column_stack((ingridx,ingridy))       
            
    
            binctr=0
            for cbin in np.unique(anglebins):
                #print "really should fix the center of rotation! to be the storm center"
                rotx=ingridx*np.cos(angs[binctr])+ingridy*np.sin(angs[binctr])
                roty=-ingridx*np.sin(angs[binctr])+ingridy*np.cos(angs[binctr])
                rotgrid=np.column_stack((rotx,roty))
                delaunay=sp.spatial.qhull.Delaunay(rotgrid)
                delarray[i].append(delaunay)
                interp=sp.interpolate.LinearNDInterpolator(delaunay,passrain.flatten(),fill_value=0.)
                tpass=np.reshape(interp(outgrid),rainprop.subdimensions)
                if rescaletype=='stochastic':
                    sys.exit("WARNING: rotation + intensity-based resampling = not tested!")
                    #whichrain[np.logical_and(whichstorms==i,anglebins==cbin)]=RainyDay.SSTalt(tpass,whichx[np.logical_and(whichstorms==i,anglebins==cbin)],whichy[np.logical_and(whichstorms==i,anglebins==cbin)],trimmask,maskheight,maskwidth,intense_data=intense_dat,durcheck=durcorrection)*rainprop.timeres/60./mnorm
                #else:
                    #whichrain[np.logical_and(whichstorms==i,anglebins==cbin)]=RainyDay.SSTalt(tpass,whichx[np.logical_and(whichstorms==i,anglebins==cbin)],whichy[np.logical_and(whichstorms==i,anglebins==cbin)],trimmask,maskheight,maskwidth,durcheck=durcorrection)*rainprop.timeres/60./mnorm     
                binctr=binctr+1
        else:
            for pt in np.arange(0,whichx.shape[3]):
                if rescaletype=='stochastic' and areatype.lower()!='pointlist' and areatype.lower!='point':                    
                    temprain,whichmultiplier[whichstorms==i,pt],whichstep=RainyDay.SSTalt(passrain,whichx[whichstorms==i,pt],whichy[whichstorms==i,pt],trimmask,maskheight,maskwidth,intensemean=intensemean,intensestd=intensestd,intensecorr=intensecorr,homemean=homemean,homestd=homestd,durcheck=durcorrection)
                    whichrain[whichstorms==i,pt]=temprain*rainprop.timeres/60./mnorm    
                elif rescaletype=='deterministic' and areatype.lower()!='pointlist' and areatype.lower()!='point':
                    temprain,whichmultiplier[whichstorms==i,pt],whichstep=RainyDay.SSTalt(passrain,whichx[whichstorms==i,pt],whichy[whichstorms==i,pt],trimmask,maskheight,maskwidth,intensemean=intensemean,homemean=homemean,durcheck=durcorrection)
                    whichrain[whichstorms==i,pt]=temprain*rainprop.timeres/60./mnorm 
                elif areatype.lower()!='pointlist' and areatype.lower()!='point' and rescaletype=='none':
                    temprain,whichstep[whichstorms==i,pt]=RainyDay.SSTalt(passrain,whichx[whichstorms==i,pt],whichy[whichstorms==i,pt],trimmask,maskheight,maskwidth,durcheck=durcorrection)                
                    whichrain[whichstorms==i,pt]=temprain*rainprop.timeres/60./mnorm 
                elif areatype.lower()=='pointlist':
                    if rescaletype=='deterministic':
                        homemeanpt=intensemean[yind_list[pt],xind_list[pt]]
                        temprain,whichmultiplier[whichstorms==i,pt],_=RainyDay.SSTalt_singlecell(passrain,whichx[whichstorms==i,pt],whichy[whichstorms==i,pt],trimmask,1,1,durcheck=durcorrection,intensemean=intensemean,homemean=homemeanpt)
                    elif rescaletype=='dimensionless':
                        homemeanpt=atlas_regridded[yind_list[pt],xind_list[pt]]
                        temprain,whichmultiplier[whichstorms==i,pt],_=RainyDay.SSTalt_singlecell(passrain,whichx[whichstorms==i,pt],whichy[whichstorms==i,pt],trimmask,1,1,durcheck=durcorrection,intensemean=atlas_regridded,homemean=homemeanpt)
                    elif rescaletype=='stochastic':
                        homemeanpt=intensemean[yind_list[pt],xind_list[pt]]
                        homestdpt=intensestd[yind_list[pt],xind_list[pt]]
                         
                        intensemeanpt=np.mean(intenserain,axis=0)
                        intensestdpt=np.std(intenserain,axis=0) 
                        
                        intensemeanpt=RainyDay.mysmoother(intensemeanpt,sigma=[smoothsig,smoothsig])
                        intensestdpt=RainyDay.mysmoother(intensestdpt,sigma=[smoothsig,smoothsig])
                         
                        intensecorrpt=np.empty((ylen_wmask,xlen_wmask),dtype='float32')    
                        intensecorrpt=RainyDay.intense_corrloop(intenserain,intensecorrpt,intenserain[:,yind_list[pt],xind_list[pt]],xlen_wmask,ylen_wmask,mnorm,domainmask)    
                        intensecorrpt=RainyDay.mysmoother(intensecorrpt,sigma=[smoothsig,smoothsig])
                         
                        intensecorrpt[np.isneginf(intensecorrpt)]=1.0
                    
                        intensemeanpt[np.isnan(intensemeanpt)]=homemeanpt
                        intensestdpt[np.isnan(intensestdpt)]=0.
                        intensecorrpt[np.isnan(intensecorrpt)]=1.0
                        temprain,whichmultiplier[whichstorms==i,pt],_=RainyDay.SSTalt_singlecell(passrain,whichx[whichstorms==i,pt],whichy[whichstorms==i,pt],trimmask,maskheight,maskwidth,intensemean=intensemeanpt,intensestd=intensestdpt,intensecorr=intensecorrpt,homemean=homemeanpt,homestd=homestdpt,durcheck=durcorrection)
                    else:
                        temprain,_=RainyDay.SSTalt_singlecell(passrain,whichx[whichstorms==i,pt],whichy[whichstorms==i,pt],trimmask,maskheight,maskwidth,durcheck=durcorrection)

                    whichrain[whichstorms==i,pt]=temprain*rainprop.timeres/60.   

                   
                elif areatype.lower()=='point':
                    
                    if rescaletype=='deterministic':
                        homemeanpt=intensemean[ymin,xmin]
                        temprain,whichmultiplier[whichstorms==i,pt],whichstep[whichstorms==i,pt]=RainyDay.SSTalt_singlecell(passrain,whichx[whichstorms==i,pt],whichy[whichstorms==i,pt],trimmask,1,1,durcheck=durcorrection,intensemean=intensemean,homemean=homemeanpt)
                        whichrain[whichstorms==i,pt]=temprain*rainprop.timeres/60.
                    elif rescaletype=='dimensionless':
                        homemeanpt=atlas_regridded[ymin,xmin]
                        temprain,whichmultiplier[whichstorms==i,pt],whichstep[whichstorms==i,pt]=RainyDay.SSTalt_singlecell(passrain,whichx[whichstorms==i,pt],whichy[whichstorms==i,pt],trimmask,1,1,durcheck=durcorrection,intensemean=atlas_regridded,homemean=homemeanpt)
                        whichrain[whichstorms==i,pt]=temprain*rainprop.timeres/60.
                    elif rescaletype=='stochastic':
                        homemeanpt=intensemean[ymin,xmin]
                        homestdpt=intensemean[ymin,xmin]
                        temprain,whichmultiplier[whichstorms==i,pt],whichstep[whichstorms==i,pt]=RainyDay.SSTalt(passrain,whichx[whichstorms==i,pt],whichy[whichstorms==i,pt],trimmask,1,1,intensemean=intensemean,intensestd=intensestd,intensecorr=intensecorr,homemean=homemeanpt,homestd=homestdpt,durcheck=durcorrection)
                        whichrain[whichstorms==i,pt]=temprain*rainprop.timeres/60. 
                    else:
                        temprain,whichstep[whichstorms==i,pt]=RainyDay.SSTalt_singlecell(passrain,whichx[whichstorms==i,pt],whichy[whichstorms==i,pt],trimmask,1,1,durcheck=durcorrection)
                        whichrain[whichstorms==i,pt]=temprain*rainprop.timeres/60. 

    if areatype.lower()=='pointlist' or areatype.lower()=='point':
        if len(arfval)==1 and np.isclose(arfval[0],1.):
            sortrain=whichrain*arfval
        else:
            if len(arfval)==1:
                arfrand=np.random.exponential(arfval[0]-1.,size=whichrain.shape)+1.
                arflimit=sp.stats.expon.ppf(0.90,1./(arfval[0]-1.))+1.
                arfmed=sp.stats.expon.ppf(0.5,1./(arfval[0]-1.))+1.
            elif len(arfval)==2:
                arfrand=np.random.gamma(shape=arfval[0],scale=arfval[1],size=whichrain.shape)
                arflimit=sp.stats.gamma.ppf(0.90,a=arfval[0],scale=arfval[1])
                arfmed=sp.stats.gamma.ppf(0.5,a=arfval[0],scale=arfval[1])
            #arfrand[arfrand>1.5]=1.0
            arfrand[arfrand>arflimit]=arfmed
            whichrain=np.multiply(whichrain,arfrand)
            
            
    # modified by DBW to account for k=0 situations, 9/22/2022
    nostorm_index=np.equal(whichstorms,-9999)
    whichrain[nostorm_index]=-9999.
    whichx[nostorm_index]=-9999
    whichy[nostorm_index]=-9999
    try:
        whichstep[nostorm_index]=-9999
    except NameError:
        pass
    try:
        whichmultiplier[nostorm_index]=-9999.
    except NameError:
        pass  
    try:
        whichtimeind[nostorm_index]=-9999
    except NameError:
        pass  
            
    
    # HERE ARE THE ANNUAL MAXIMA!!!
    if areatype.lower()=="pointlist":
        sortrain=np.empty((whichrain.shape[1],whichrain.shape[2],npoints_list),dtype='float32')
        #if rescaletype=='stochastic' or rescaletype=='deterministic':
        #    sortmultiplier=np.empty((whichrain.shape[1],whichrain.shape[2],npoints_list),dtype='float32')
        for pt in np.arange(0,whichx.shape[3]):
            
            if calctype.lower()=='ams':
                maxrain=np.nanmax(whichrain[:,:,:,pt],axis=0)
                
            elif calctype.lower()=='pds':
                temprain=np.squeeze(whichrain[:,:,:,pt])
                temppds=temprain.reshape(-1, temprain.shape[-1])
                maxrain=np.sort(temppds,axis=0)[-nsimulations:,:]
            
            
            # PULL OUT THE CORRESPONDING TRANSPOSITION INFORMATION
            maxind=np.nanargmax(whichrain[:,:,:,pt],axis=0)
            
            #added for k=0 catching, DBW 9/22/22
            maxind[np.equal(maxrain,-9999.)]=-9999
        
            
            # THIS ISN'T VERY ELEGANT
            maxx=np.empty((maxind.shape),dtype="int32")
            maxy=np.empty((maxind.shape),dtype="int32")
            maxstorm=np.empty((maxind.shape),dtype="int32")
            maxx[:]=-9999
            maxy[:]=-9999
            maxstorm[:]=-9999
            #if rescaletype=='stochastic' or rescaletype=='deterministic':
            #    maxmultiplier=np.empty((maxind.shape),dtype="float32") 
                
            for i in range(0,np.max(ncounts)):
                maxx[maxind==i]=whichx[i,maxind==i,pt]
                maxy[maxind==i]=whichy[i,maxind==i,pt]
                maxstorm[maxind==i]=whichstorms[i,maxind==i]
             #   if rescaletype=='stochastic' or rescaletype=='deterministic':
             #       maxmultiplier[maxind==i]=whichmultiplier[i,maxind==i,pt]
            
            
            # RANK THE STORMS BY INTENSITY AND ASSIGN RETURN PERIODS
            exceedp=np.linspace(1,1./nsimulations,nsimulations)
            
            returnperiod=1/exceedp
            sortind=np.argsort(maxrain,axis=0)
            sortrain[:,:,pt]=np.sort(maxrain,axis=0)
            #if rescaletype=='stochastic' or rescaletype=='deterministic':
            #    sortmultiplier[:,:,pt]=maxmultiplier[sortind]
                
                
        if alllevels==False:
            reducedlevind=[]
            for i in range(0,len(speclevels)):
                reducedlevind.append(RainyDay.find_nearest(returnperiod,speclevels[i]))  
            
            returnperiod=returnperiod[reducedlevind]
            sortrain=sortrain[reducedlevind,:]
            exceedp=exceedp[reducedlevind]
                
    else:
        if calctype=='ams' or calctype=='pds':          # this isn't very elegant!
            if calctype.lower()=='ams':
                maxrain=np.nanmax(whichrain[:,:,:,pt],axis=0)
                maxind=np.nanargmax(whichrain,axis=0)
            elif calctype.lower()=='pds':
                temprain=np.squeeze(whichrain)
                temppds=temprain.reshape(-1, temprain.shape[-1])
                maxrain=np.sort(temppds,axis=0)[-nsimulations:,:]
                maxind=np.nanargmax(whichrain,axis=0)
            
            #added for k=0 catching, DBW 9/22/22
            maxind[np.equal(maxrain,-9999.)]=-9999
            
            #elif calctype.lower()=='npyear':
            #    temprain=np.squeeze(whichrain)
            #    maxrain=np.sort(temprain,axis=0)
            #    maxind=np.argsort(temprain,axis=0)
                
                #maxrain=np.sort(temprain,axis=0)[-5:,:]
                #maxind=np.argsort(temprain,axis=0)[-5:,:]
      
            
            # HERE THE OPTIONAL USER SPECIFIED INTENSITY DISTRIBUTION IS APPLIED    
            if userdistr.all()!=False:
                rvs=sp.stats.genextreme.rvs(userdistr[2],loc=userdistr[0],scale=userdistr[1],size=maxrain.shape).astype('float32')
                maxrain=maxrain*rvs
                maxrain[np.equal(maxind,-9999)]=-9999.
                
            # PULL OUT THE CORRESPONDING TRANSPOSITION INFORMATION
            
            
            # THIS ISN'T VERY ELEGANT
            maxx=np.empty((maxind.shape),dtype="int32")
            maxy=np.empty((maxind.shape),dtype="int32")
            maxx[:]=-9999
            maxy[:]=-9999
            maxstorm=np.empty((maxind.shape),dtype="int32")
            maxstorm[:]=-9999
            if arfcorrection:
                maxstep=np.empty((maxind.shape),dtype="int32")
                maxstep[:]=-9999
            if rotation:
                maxangles=np.empty((maxind.shape),dtype="float32")
                sortangle=np.empty((maxind.shape),dtype="float32")
                maxangles[:]=-9999.
                sortangle[:]=-9999.
            if rescaletype=='stochastic' or rescaletype=='deterministic' or rescaletype=='dimensionless':
                maxmultiplier=np.empty((maxind.shape),dtype="float32") 
                sortmultiplier=np.empty((maxind.shape),dtype="float32")
                maxmultiplier[:]=-9999.
                sortmultiplier[:]=-9999.
                
            for i in range(0,np.max(ncounts)):
                maxx[maxind==i]=np.squeeze(whichx[i,np.squeeze(maxind==i)])
                maxy[maxind==i]=np.squeeze(whichy[i,np.squeeze(maxind==i)])
                maxstorm[maxind==i]=np.squeeze(whichstorms[i,np.squeeze(maxind==i)])
                if arfcorrection:
                    maxstep[maxind==i]=np.squeeze(whichstep[i,np.squeeze(maxind==i)])
                
                if rotation:
                    maxangles[maxind==i]=randangle[i,maxind==i]
                if rescaletype=='stochastic' or rescaletype=='deterministic' or rescaletype=='dimensionless':
                    maxmultiplier[maxind==i]=np.squeeze(whichmultiplier[i,maxind==i])
    #            elif calctype.lower()=='npyear':
    #                sys.exit("having problems here")
    #                for stm in range(0,nperyear):
    #                    maxx[maxind==i]=np.squeeze(whichx[i,np.squeeze(maxind==i)])
    #                    maxy[stm,maxind==i]=np.squeeze(whichy[i,np.squeeze(maxind==i)])
    #                    maxstorm[stm,maxind==i]=np.squeeze(whichstorms[i,np.squeeze(maxind==i)])
    #                    if rotation:
    #                        print("Warning: We haven't tested ROTATION with NPERYEAR!")
    #                        maxangles[stm,maxind==i]=randangle[stm,i,maxind==i]
    #                    if rescaletype=='stochastic' or rescaletype=='deterministic' or rescaletype=='dimensionless':
    #                        maxmultiplier[stm,maxind==i]=np.squeeze(whichmultiplier[stm,i,maxind==i])
                    
                
            
            # RANK THE STORMS BY INTENSITY AND ASSIGN RETURN PERIODS
            exceedp=np.linspace(1,1./nsimulations,nsimulations)
            returnperiod=1/exceedp
            #rp_pds=1./(1.-np.exp(-1./(returnperiod)))
            #rp_pds=1./np.log(returnperiod/(returnperiod-1))
            sortind=np.argsort(maxrain,axis=0)
            sortrain=np.sort(maxrain,axis=0)
            sortx=np.empty((maxind.shape),dtype="int32")
            sorty=np.empty((maxind.shape),dtype="int32")
            sortstorms=np.empty((maxind.shape),dtype="int32")
            if arfcorrection:
                sortstep=np.empty((maxind.shape),dtype="int32")
            
            
            for i in range(0,nrealizations):
                sortx[:,i]=maxx[sortind[:,i],i]
                sorty[:,i]=maxy[sortind[:,i],i]
                sortstorms[:,i]=maxstorm[sortind[:,i],i]
                if arfcorrection:
                    sortstep[:,i]=maxstep[sortind[:,i],i]
                if rotation:
                    sortangle[:,i]=maxangles[sortind[:,i],i]
                if rescaletype=='stochastic' or rescaletype=='deterministic' or rescaletype=='dimensionless':
                    sortmultiplier[:,i]=maxmultiplier[sortind[:,i],i]
            
                
            # FIND THE TIMES:
            if arfcorrection==False:
                sorttimes=np.zeros((maxind.shape[0],maxind.shape[1],cattime.shape[1]),dtype="datetime64[m]")
            else:       # using ARFANALYSIS
                sorttimes=np.zeros((maxind.shape[0],maxind.shape[1],int(duration*60/rainprop.timeres)),dtype="datetime64[m]")
            whichorigstorm=np.zeros((maxind.shape[0],maxind.shape[1]),dtype='int32')
            for i in range(0,nstorms):
                if np.sum(np.squeeze(sortstorms==i))>0:
                    if arfcorrection==False:
                        sorttimes[np.squeeze(sortstorms==i),:]=cattime[i,:]
                        
                    else:  # using ARFANALYSIS
                        tstep=sortstep[np.squeeze(sortstorms==i),:]
                        if len(tstep)>1:
                            tempstep=np.squeeze(sortstep[np.squeeze(sortstorms==i),:])
                            temptime=np.empty((tempstep.shape[0],int(duration*60/rainprop.timeres)),dtype="datetime64[m]")
                            for j in range(0,tempstep.shape[0]):
                                temptime[j,:]=cattime[i,tempstep[j]:tempstep[j]+int(duration*60/rainprop.timeres)]
                        else:
                            tempstep=sortstep[np.squeeze(sortstorms==i),:][0]
                            temptime=np.empty((tempstep.shape[0],int(duration*60/rainprop.timeres)),dtype="datetime64[m]")
                            for j in range(0,tempstep.shape[0]):
                                temptime[j,:]=cattime[i,tempstep[j]:tempstep[j]+int(duration*60/rainprop.timeres)]
    
                        sorttimes[np.squeeze(sortstorms==i),:]=temptime
                    whichorigstorm[np.squeeze(sortstorms==i)]=modstormsno[i]+1
                else:
                    continue
                
                
            if alllevels==False:
                reducedlevind=[]
                for i in range(0,len(speclevels)):
                    reducedlevind.append(RainyDay.find_nearest(returnperiod,speclevels[i]))  
                
                returnperiod=returnperiod[reducedlevind]
                sortrain=sortrain[reducedlevind,:]
                sortstorms=sortstorms[reducedlevind,:]
                sorttimes=sorttimes[reducedlevind,:]
                exceedp=exceedp[reducedlevind]
                sortx=sortx[reducedlevind,:]
                sorty=sorty[reducedlevind,:]
                if arfcorrection:
                    sortstep=sortstep[reducedlevind,:]

                whichorigstorm=whichorigstorm[reducedlevind,:]
                if rotation:    
                    sortangle=sortangle[reducedlevind,:]
                if rescaletype=='stochastic' or rescaletype=='deterministic' or rescaletype=='dimensionless':        
                    sortmultiplier=sortmultiplier[reducedlevind,:]
            
        
        #################################################################################
        # STEP 2a (OPTIONAL): Find the single storm maximized storm rainfall-added DBW 7/19/2017
        #################################################################################    
        
        nanmask=deepcopy(trimmask)
        nanmask[np.isclose(nanmask,0.)]=np.nan
        nanmask[np.isclose(nanmask,0.)==False]=1.0
        
        if deterministic:
            print("finding maximizing precipitation...")
            
            max_trnsx=catx[-1]
            max_trnsy=caty[-1]
            if rotation==False:
                # there is some small bug that I don't understand either here or in the storm catalog creation, in which maxstm_avgrain will not exactly match catmax[-1] unless areatype is a point
                maxstm_rain=np.multiply(catrain[-1,:,max_trnsy:(max_trnsy+maskheight),max_trnsx:(max_trnsx+maskwidth)],nanmask)
                maxstm_avgrain=np.nansum(np.multiply(catrain[-1,:,max_trnsy:(max_trnsy+maskheight),max_trnsx:(max_trnsx+maskwidth)],trimmask))/mnorm
                maxstm_ts=np.nansum(np.multiply(maxstm_rain,trimmask)/mnorm,axis=(1,2))
                maxstm_time=cattime[-1,:]
            else:  
                prevmxstm=0.
                maxstm_rain=np.empty((catrain.shape[1],nanmask.shape[0],nanmask.shape[1]),dtype='float32')
                for i in range(0,nstorms):
                    passrain=np.nansum(catrain[i,:],axis=0)
                    xctr=catx[i]+maskwidth/2.
                    yctr=caty[i]+maskheight/2.
                    xlinsp=np.linspace(-xctr,rainprop.subdimensions[1]-xctr,rainprop.subdimensions[1])
                    ylinsp=np.linspace(-yctr,rainprop.subdimensions[0]-yctr,rainprop.subdimensions[0])
                    ingridx,ingridy=np.meshgrid(xlinsp,ylinsp)
                    ingridx=ingridx.flatten()
                    ingridy=ingridy.flatten()
                    outgrid=np.column_stack((ingridx,ingridy))       
                
                    for tempang in angbins:
                        #print "really should fix the center of rotation! to be the storm center"
                        rotx=ingridx*np.cos(tempang)+ingridy*np.sin(tempang)
                        roty=-ingridx*np.sin(tempang)+ingridy*np.cos(tempang)
                        rotgrid=np.column_stack((rotx,roty))
                        delaunay=sp.spatial.qhull.Delaunay(rotgrid)
                        interp=sp.interpolate.LinearNDInterpolator(delaunay,passrain.flatten(),fill_value=0.)
                        train=np.reshape(interp(outgrid),rainprop.subdimensions)
                        temp_maxstm_avgrain=np.nansum(np.multiply(train[max_trnsy:(max_trnsy+maskheight),max_trnsx:(max_trnsx+maskwidth)],trimmask))/mnorm
                        if temp_maxstm_avgrain>prevmxstm:
                            maxstm_avgrain=temp_maxstm_avgrain
                            prevmxstm=maxstm_avgrain
                            maxstm_time=cattime[-i,:]
                            
                            for k in range(0,len(maxstm_time)):
                                interp=sp.interpolate.LinearNDInterpolator(delaunay,catrain[i,k,:].flatten(),fill_value=0.)
                                maxstm_rain[k,:]=np.reshape(interp(outgrid),rainprop.subdimensions)[max_trnsy:(max_trnsy+maskheight),max_trnsx:(max_trnsx+maskwidth)]
                            maxstm_rain=np.multiply(maxstm_rain,nanmask)
                            maxstm_ts=np.nansum(np.multiply(maxstm_rain,trimmask)/mnorm,axis=(1,2))
        	
        
        
    #################################################################################
    # STEP 3 (OPTIONAL): RAINFALL FREQUENCY ANALYSIS
    #################################################################################
    
    if calctype!='npyear':
        print("preparing frequency analysis...")
    
        if areatype.lower()=='pointlist':
            spreadmean=np.nanmean(sortrain,1)
            
            if spreadtype=='ensemble':
                spreadmin=np.nanmin(sortrain,axis=1)
                spreadmax=np.nanmax(sortrain,axis=1)    
            else:
                spreadmin=np.percentile(sortrain,(100-quantilecalc)/2,axis=1)
                spreadmax=np.percentile(sortrain,quantilecalc+(100-quantilecalc)/2,axis=1)
            
            if spreadmean.shape[0]!=returnperiod.shape[0] or spreadmax.shape[0]!=returnperiod.shape[0] or spreadmean.shape[1]!=ptlatlist.shape[0]:
                sys.exit("There is some dimension inconsistency in the pointlist scheme!")
            
            #fmean=open(FreqFile_mean,'w')
            #fmin=open(FreqFile_min,'w')
            #fmax=open(FreqFile_max,'w')
            
            #fmean.write('#prob.exceed,returnperiod,meanrain\n'+ptlistname+'\n')
            #fmax.write('#prob.exceed,returnperiod,maxrain\n'+ptlistname+'\n')
            #fmin.write('#prob.exceed,returnperiod,minrain\n'+ptlistname+'\n')
    
            freqanalysis_mean=np.column_stack((exceedp,returnperiod,spreadmean))
            freqanalysis_min=np.column_stack((exceedp,returnperiod,spreadmin))
            freqanalysis_max=np.column_stack((exceedp,returnperiod,spreadmax))
            
            outlat_line=np.append([-999.,-999.],ptlatlist)
            outlon_line=np.append([-999.,-999.],ptlonlist)
            coordline=np.row_stack((outlat_line,outlon_line))
            freqanalysis_mean=np.row_stack((coordline,freqanalysis_mean))
            freqanalysis_min=np.row_stack((coordline,freqanalysis_min))
            freqanalysis_max=np.row_stack((coordline,freqanalysis_max))
            
            np.savetxt(FreqFile_mean,freqanalysis_mean,delimiter=',',header='prob.exceed,returnperiod,meanrain',fmt='%6.2f',comments='#',footer=ptlistname)
            np.savetxt(FreqFile_min,freqanalysis_min,delimiter=',',header='prob.exceed,returnperiod,minrain',fmt='%6.2f',comments='#',footer=ptlistname)
            np.savetxt(FreqFile_max,freqanalysis_max,delimiter=',',header='prob.exceed,returnperiod,maxrain',fmt='%6.2f',comments='#',footer=ptlistname)
            
        else:
            if spreadtype=='ensemble':
                spreadmin=np.nanmin(sortrain,axis=1)
                spreadmax=np.nanmax(sortrain,axis=1)    
            else:
                spreadmin=np.percentile(sortrain,(100-quantilecalc)/2,axis=1)
                spreadmax=np.percentile(sortrain,quantilecalc+(100-quantilecalc)/2,axis=1)
        
            freqanalysis=np.column_stack((exceedp,returnperiod,spreadmin,np.nanmean(sortrain,1),spreadmax))
            
            np.savetxt(FreqFile,freqanalysis,delimiter=',',header='prob.exceed,returnperiod,minrain,meanrain,maxrain',fmt='%6.2f',comments='')
            
            import matplotlib.patches as mpatches
            from matplotlib.font_manager import FontProperties
            from matplotlib import pyplot as plt
           # warnings.filterwarnings('ignore')
            fontP = FontProperties()
            fontP.set_size('xx-small')
            fig, ax = plt.subplots(1)
            line1, = plt.plot(exceedp[exceedp<=0.5], RainyDay.np.nanmean(sortrain,1)[exceedp<=0.5], lw=1, label='Average', color='blue')
        
            ax.fill_between(exceedp[exceedp<=0.5], spreadmin[exceedp<=0.5], spreadmax[exceedp<=0.5], facecolor='dodgerblue', alpha=0.5,label='Ensemble Variability')
            blue_patch = mpatches.Patch(color='dodgerblue', label='Spread')
            plt.legend(handles=[line1,blue_patch],loc='lower right',prop = fontP)
        
            if np.nanmax(spreadmax[exceedp<=0.5])<10.:
                upperlimit=10.
            elif np.nanmax(spreadmax[exceedp<=0.5])<100.:
                upperlimit=100.
            elif np.nanmax(spreadmax[exceedp<=0.5])<1000.:
                upperlimit=1000.
            else:
                upperlimit=10000.
                
            if np.nanmin(spreadmin[exceedp<=0.5])<1.:
                lowerlimit=0.1
            elif np.nanmin(spreadmin[exceedp<=0.5])<10.:
                lowerlimit=1 
            elif np.nanmin(spreadmin[exceedp<=0.5])<100.:
                lowerlimit=10. 
            else:
                lowerlimit=100.
                    
                    
            plt.ylim(lowerlimit,upperlimit)
            ax.set_xlabel('Annual Exceed. Prob. [-]\n1/(Return Period) [year]')
            ax.set_ylabel('Precip. Depth [mm]')
            ax.set_yscale('log')
            ax.set_xscale('log')
            plt.gca().invert_xaxis()
            ax.grid()
            plt.tight_layout()
            plt.savefig(fullpath+'/'+scenarioname+'_FrequencyAnalysis.png',dpi=250)
            plt.close('all')
    else:
        print("You selected NPERYEAR greater than 1. RainyDay isn't configured to do a rainfall frequency analysis for this, but will write output scenarios. If you want a rainfall frequency analysis, set NPERYEAR to '1','false' or omit it from the .sst file!")
       
           