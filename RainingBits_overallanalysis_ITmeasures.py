#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 10:31:22 2018

Updated May 2020

This program computes multidimensional lagged MI between lagged and current rainfall
histories, for 70 year period, for four seasons
*results in list of information measures described in Goodwell 2020 JHM, for each grid cell, year, and season

Rainfall from CPC gage-based gridded data, converted to binary according to 
%-ile threshold (after values below thresold have been removed)

consider lagged mutual information between neighboring grid cells (8 directions)
resolve vector of information values to obtain "dominant direction"

Inputs: 
precipitation files from CPC Unified gage-based dataset (annual .nc files)*
results from threshold finding code: 'Spatial_thresholds_mm.pickle' file
*note: need to change precip_folder variable to reflect location of precip data

Outputs:
pickle file called 'spatial_results_overall.pickle'
with embedded lists for rows of latitude

Notes:
*this code runs compute_info_measures function in parallel for each row of latitude
must choose number of processors in pool
*next run code CPC_overall_analyzeresults to analyze results for maps in Goodwell2020

@author: allisongoodwell

@author: allisongoodwell
"""

from netCDF4 import Dataset
import numpy as np
import pickle

from multiprocessing import Pool
import time
import gc


precip_folder = 'C:\\Users\\goodwela\\Dropbox\\UCDenver\\research\\rainfall_research\\CPC_raingage_gridded_dailydata\\original_precip_data\\'
nProcessors = 14

np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

maxdim = 9 #9 lagged neighbors (including self) and 1 current target at center

lags = [1]
N = 2 #number of bins (2 for binary pdf)
years = range(1948,2018)
nyrs = len(years)


#function to compute statistical significance, shuffled surrogates
def compute_stat_sig(Tuple):
    
    #shuffle all dimensions except last one (the target)
    ntests=10
    sig_level = 4

    Itot_vector = np.zeros(ntests)
    for j in range(0,ntests):
        
        np.random.shuffle(Tuple[:,0:-1])
                    
        pdf,edges = np.histogramdd(Tuple,bins=2)
        pdf = pdf/np.sum(pdf)
        source_dim = tuple(range(0,np.shape(Tuple)[1]-1))
        tar_dim = np.shape(Tuple)[1]-1;
        
        m_sources = np.sum(pdf,axis=tar_dim) #should be dim-1 dimension pdf
        m_target = np.sum(pdf,axis=source_dim) #should be 1D pdf
        #pad m_sources with 1 for last index
        m_sources = np.expand_dims(m_sources, axis=np.size(np.shape(m_sources)))
             
        p_independent = m_sources * m_target
        f_ijk = np.divide(pdf,p_independent,out=np.zeros_like(pdf),where=p_independent!=0)               
        logf=np.log2(f_ijk,out=np.zeros_like(f_ijk),where=f_ijk!=0)      
        plogf = pdf*logf
        
        Itot_vector[j] = np.nansum(plogf)
    
    I_mean=np.average(Itot_vector)
    I_std=np.std(Itot_vector)  
    Icrit = I_mean + sig_level*I_std
    
    return Icrit

#function to compute seasonal pdfs and info measures
def compute_info_measures(input_list):
 
    ppt_list = input_list[0]
    years = input_list[1]
    threshvals = input_list[2]
    days = input_list[3]
    
    info_all=[]
    for s in range(0,4):
        
        #for each latitude and longitude pair, find 8 nearest neighbors                 
        ppt_list_mini = [ppt_list[y_ind][days[s],:,:] for y_ind,y in enumerate(years)]
        PPTvect = np.concatenate(ppt_list_mini)
        
        cutoff_starts = [y*len(days[s])-1 for y in range(1,len(years))]
        
        if np.isnan(PPTvect[0,1,1])==True:
            #print('skipping!')
            continue
        #PPTvect is all precip for a given coordinate, and 8 surrounding coordinates
        # size time length x latitudes x longitudes
        
        PPTvect[np.isnan(PPTvect)]=0
        
        #get binary vector of all rainfall values
        pptvect_binary = np.zeros(np.shape(PPTvect))
        
        pptvect_binary[PPTvect>threshvals[s]]=1         
    
        Xtar = pptvect_binary[1:,1,1]        
        lagged_vect = pptvect_binary[:-1,:,:]
        
        #remove season breaks from Xtar and Xlag (need to find indices)
        #lagged_vect = [lagged_vect[i,:,:] for i in range(0,len(Xtar)) if i not in cutoff_starts]
        #Xtar = [Xtar[i] for i in range(0,len(Xtar)) if i not in cutoff_starts]

        
        lagged_vect = np.delete(lagged_vect,cutoff_starts,axis=0)
        Xtar = np.delete(Xtar,cutoff_starts,axis=0)
       
        #compute 2d mutual lagged mutual information in each direction
        #compute total lagged mutual information from top 5 directions
        
        #Tuple = np.vstack((lagged_vect[:,1,2],lagged_vect[:,2,2],lagged_vect[:,2,1],lagged_vect[:,2,0],
                        #   lagged_vect[:,1,0],lagged_vect[:,0,0],lagged_vect[:,0,1],lagged_vect[:,0,2],
                        #   lagged_vect[:,1,1],Xtar))
                        
            
        
        Tuple = np.vstack((lagged_vect[:,1,2],lagged_vect[:,2,2],lagged_vect[:,2,1],lagged_vect[:,2,0],
                   lagged_vect[:,1,0],lagged_vect[:,0,0],lagged_vect[:,0,1],lagged_vect[:,0,2],
                   lagged_vect[:,1,1],Xtar))
              
        Tuple = np.transpose(Tuple)
                                
        #compute 2d MI between each entity (all neighbors and lagged self)
        #also compute joint MI between lagged self and each neighbor to current target
        I_2D=[]
        I_3D=[]
        Eta =[]
        Eta_bits=[]
        for d in range(0,9):
            
            if np.sum(Tuple[:,d])<1:
                I_2D.append(0)
                continue
            
            smallTuple = Tuple[:,[d,-1]]
            pdf,edges = np.histogramdd(smallTuple,bins=2)
            pdf = pdf/np.sum(pdf)
            
            source_dim = 0;
            tar_dim = 1;
            
            m_sources = np.sum(pdf,axis=tar_dim) #should be dim-1 dimension pdf
            m_target = np.sum(pdf,axis=source_dim) #should be 1D pdf
            
            #Hx = np.sum(m_target*np.log2(m_target))*-1
           
            #pad m_sources with 1 for last index
            m_sources = np.expand_dims(m_sources, axis=np.size(np.shape(m_sources)))
                 
            p_independent = m_sources * m_target
            f_ijk = np.divide(pdf,p_independent,out=np.zeros_like(pdf),where=p_independent!=0)               
            logf=np.log2(f_ijk,out=np.zeros_like(f_ijk),where=f_ijk!=0)     
            plogf = pdf*logf
            
            Ival = np.nansum(plogf)
            
            #check for statistical significance of Ival
            Icrit = compute_stat_sig(smallTuple)
            
            #print(Ival,Icrit)
            
            if Icrit > Ival:
                Ival=0
            
            I_2D.append(Ival)
            
            
            
        for d in range(0,9):
            
            if np.sum(Tuple[:,d])<1:
                I_3D.append(0)
                Eta_bits.append(0)
                Eta.append(0)
                continue      
            
            Tuple3d = Tuple[:,[d,-2,-1]]
            pdf,edges = np.histogramdd(Tuple3d,bins=2)
            pdf = pdf/np.sum(pdf)
            
            source_dim = tuple(range(0,2))
            tar_dim = 2;
            
            m_sources = np.sum(pdf,axis=tar_dim) #should be dim-1 dimension pdf
            m_target = np.sum(pdf,axis=source_dim) #should be 1D pdf
            
            #pad m_sources with 1 for last index
            m_sources = np.expand_dims(m_sources, axis=np.size(np.shape(m_sources)))
                 
            p_independent = m_sources * m_target
            f_ijk = np.divide(pdf,p_independent,out=np.zeros_like(pdf),where=p_independent!=0)               
            logf=np.log2(f_ijk,out=np.zeros_like(f_ijk),where=f_ijk!=0)     
            plogf = pdf*logf
            
            Ival = np.nansum(plogf)
            
            #now break this down into U_here, U_there, R and S components
            I_here = I_2D[-1] #lagged info between center cell and itself is last index of I_2D
            #I_there = I_2D[d] #ind lagged info for neighbor and current center cell
            
            #check for statistical significance of Ival (but only if I_here is 0)
            if I_here ==0:
                Icrit = compute_stat_sig(Tuple3d)
                if Icrit > Ival:
                    Ival=0
            
            I_3D.append(Ival)
            
            #determine eta = U_there + S / Itot as the information proportion above that from center itself
            Eta_proportion = np.divide(Ival-I_here,Ival,out=np.zeros_like(Ival),where=Ival!=0)
            Eta.append(Eta_proportion)
            
            Eta_bits.append(Ival-I_here)
            
        #next item: use Eta_bits to pick the top 4 U and S neighbors
        #compute total information from those 4 neighbors
        top_idx = np.argsort(Eta_bits)[-4:]
        top_index = np.append(top_idx,9)
        
        newTuple = Tuple[:,top_index]
        pdf,edges = np.histogramdd(newTuple,bins=N)
        pdf = pdf/np.sum(pdf)
               
        source_dim = tuple(range(0,4))        
        tar_dim = 4 #target dim is always the max
            
        m_sources = np.sum(pdf,axis=tar_dim) #should be dim-1 dimension pdf
        m_target = np.sum(pdf,axis=source_dim) #should be 1D pdf
            
        #pad m_sources with 1 for last index
        m_sources = np.expand_dims(m_sources, axis=np.size(np.shape(m_sources)))
                 
        p_independent = m_sources * m_target
        f_ijk = np.divide(pdf,p_independent,out=np.zeros_like(pdf),where=p_independent!=0)               
        logf=np.log2(f_ijk,out=np.zeros_like(f_ijk),where=f_ijk!=0)     
        plogf = pdf*logf
            
        I4=np.nansum(plogf) #I4 is total info from 4 sources
        
        #angles = np.asfarray([0, 45, 90, 135, 180, 225, 270, 315])
        #last item: opposing joint MI measures (N-S, E-W, NE-SW, NW-SE)
        pairs = [[0,4],[1,5],[2,6],[3,7]] #E-W, NE-SW, N-S, NW-SE pairs of indices
        
        U1=[]
        U2=[]
        R=[]
        S=[]
        Idirs=[]
        
        for p in pairs:
            p.append(9)
    
            newTuple = Tuple[:,p]
            pdf,edges = np.histogramdd(newTuple,bins=N)
            pdf = pdf/np.sum(pdf)
            
            
            source_dim = tuple(range(0,2))        
            tar_dim = 2 #target dim is always the max
                
            m_sources = np.sum(pdf,axis=tar_dim) #should be dim-1 dimension pdf
            m_target = np.sum(pdf,axis=source_dim) #should be 1D pdf
                
            #pad m_sources with 1 for last index
            m_sources = np.expand_dims(m_sources, axis=np.size(np.shape(m_sources)))
                     
            p_independent = m_sources * m_target
            f_ijk = np.divide(pdf,p_independent,out=np.zeros_like(pdf),where=p_independent!=0)               
            logf=np.log2(f_ijk,out=np.zeros_like(f_ijk),where=f_ijk!=0)        
            plogf = pdf*logf
                
            Itotal = np.nansum(plogf)
            
            #check for statistical significance of Itotal (but only if I_2d components =0)
            if I_2D[p[0]]==0 and I_2D[p[1]]==0:
                Icrit = compute_stat_sig(newTuple)
                if Icrit>Itotal:
                    Itotal =0
            
            Idirs.append(Itotal)
                 
            #also want to break down into U,R,S components from opposite directions
            H = -1 * m_target[0] * np.log2(m_target[0]) - 1 * m_target[1]*np.log2(m_target[1])
    
            Is = Itotal/H
        
            R_MMI= np.min(np.stack((I_2D[p[0]],I_2D[p[1]])))
            R_min = np.max((-Itotal + I_2D[p[0]] + I_2D[p[1]],0.00))
    
            rval = np.max(R_min + Is* (R_MMI-R_min),0)
    
            R.append(rval)            
            U1.append(I_2D[p[0]] - rval)
            U2.append(I_2D[p[1]] - rval)
            S.append(Itotal - I_2D[p[0]]-I_2D[p[1]]+rval)
                                      
        infodict = {'I_3d':I_3D,'I_2d':I_2D,'Eta':Eta_bits,'I_top4':I4, 'R':R, 'S':S, 'UN':U1,'US':U2,'I_directions':Idirs}
        
        
        info_all.append(infodict) #list of dictionaries (one for each season)
     
    return info_all

#%%
if __name__ == '__main__':
    __spec__=None
    ppt_list=[]
    for y in years:
        
        string = precip_folder + 'precip.V1.0.'+str(y)+'.nc'
        dataset = Dataset(string)
        
        pptdata = np.asarray(dataset.variables['precip'][:])
        lat = dataset.variables['lat'][:]
        lon = dataset.variables['lon'][:]
        t = dataset.variables['time'][:]
            
        pptdata[pptdata<0]=np.nan
        
        ppt_list.append(pptdata)
        
    #determine size of a single string of ppt data
    ppt_list_mini = [ppt_list[x][:,0,0] for x in range(0,np.size(ppt_list))]
    PPTvect = np.concatenate(ppt_list_mini)
    
    npts = len(PPTvect)-1
    
    
    ndays_per_month = np.asfarray([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    ndays_cum = np.cumsum(ndays_per_month)
    
    days=[]
    days.append([i for i in range(0,365) if i<ndays_cum[1] or i >=ndays_cum[10]])
    days.append([i for i in range(0,365) if i>=ndays_cum[1] and i < ndays_cum[4]])
    days.append([i for i in range(0,365) if i>=ndays_cum[4] and i < ndays_cum[7]])
    days.append([i for i in range(0,365) if i>=ndays_cum[7] and i < ndays_cum[10]])
    
    del days[0][-1]
        
    print('data loaded....')
    
    #%% determine treshold rainfall value for each season and grid cell
    
    thresh_file='Spatial_thresholds_mm.pickle'                  
    f_myfile = open(thresh_file,'rb')
    threshvals = pickle.load(f_myfile)
    f_myfile.close()        
                
    #%% 
    
    maxinfoindex=[]
    AllResults=[]
    Long_Lists=[]
    Lat=[]
    for la_ind,la in enumerate(lat):
        
        #if np.mod(la_ind,10)==0:
        print('latitude = '+str(la))
        
        if la_ind<1 or la == np.max(lat):
            continue

        data_list=[]
        lon_list=[]
        lat_list=[] 
        item_list=[]    
        for lo_ind,lo in enumerate(lon):
            
            #t1 = time.clock()
            
            ppt_test = ppt_list[0][0,la_ind,lo_ind]
            
            if lo_ind<1 or lo==np.max(lon) or np.isnan(ppt_test)==True:
                continue
            
            print('longitude = '+str(lo) + ' latitude = ' +str(la))
            
            
            ppt_data = [ppt_list[y_ind][:,(la_ind-1):(la_ind+2),(lo_ind-1):(lo_ind+2)] for y_ind,y in enumerate(years)]
                   
            item_list = [ppt_data, years, threshvals[la_ind,lo_ind,:],days]
           
            data_list.append(item_list)
            lon_list.append(lo)
            #now want to move all the below stuff into a function that takes an element from data_list
    
        gc.collect()    
        t1 = time.clock()
        print('going into metrics with '+ str(len(lon_list))+ ' for latitude '+ str(la))    
        p = Pool(nProcessors)
        results = p.map(compute_info_measures,data_list) #
        p.close
        p.join
      
        t2 = time.clock()  
        print('time elapsed for latitude row = '+ str(t2-t1))    
        
        gc.collect()
        
#        for r_ind,r in enumerate(results):
#            print('index = '+ str(r_ind))
#            print('longitude ='+str(lon_list[r_ind]))
#            print('I value = '+ str(r[0]['I_top4']))
                        
        AllResults.append(results)   
        Long_Lists.append(lon_list)
        Lat.append(la)
                
    site_filenamesave = 'spatial_results_overall.pickle'
    f_myfile = open(site_filenamesave,'wb')
    pickle.dump(AllResults,f_myfile)
    pickle.dump(Long_Lists,f_myfile)
    pickle.dump(Lat,f_myfile)
    f_myfile.close()   
    
    print('saving results pickle file')
        
    #now need to do some filling to obtain maps from lists of dictionaries with results
    #%%
    #each item in AllResults is associated with a latitude from Lat
    #then each item in AllResults[0] is associated with a longitude from  Long_lists[0]
    #each item in AllResults[0][0] is a season 
    #AllResults[0][0][0]['I_2d'] should access a specific measure
#    lats = np.asfarray(lat)
#    lons = np.asfarray(lon)
#    
#    for a_ind,A in enumerate(AllResults): #A is list of longitudes
#        for b_ind,B in enumerate(A): #B is single longitude
#            for s, info in enumerate(B): #info is results structure for single season s
#                
#                #need to figure out what these are:
#                lat_value = Lat[a_ind]
#                lon_value = Long_Lists[a_ind][b_ind]
#                
#                la_ind = np.where(lats==lat_value)
#                lo_ind = np.where(lons==lon_value)
#      
#                I_2dvals[la_ind,lo_ind,s,:]=info['I_2d']
#                I_3dvals[la_ind,lo_ind,s,:]=info['I_3d']
#                Eta_vals[la_ind,lo_ind,s,:]=info['Eta']
#                
#                #compute primary direction (angle) and strength of info to center based on 2d metrics
#                Ivals = info['I_2d']
#                Ivals = Ivals[0:8] #index 8 is lagged middle cell itself, 9 is current value of middle cell
#                
#                I_xval = np.sum([np.cos(angles[i])*Ivals[i] for i in range(0,8)])
#                I_yval = np.sum([np.sin(angles[i])*Ivals[i] for i in range(0,8)])
#                
#                I_vectorstrength[la_ind,lo_ind,s]=np.sqrt(I_xval**2 + I_yval**2)
#                I_mainangle[la_ind,lo_ind,s] = np.arctan2(I_yval,I_xval)
#                
#                #compute primary direction (angle) and strength of info to center based on TE metrics
#                Etavals = info['Eta']
#                Etavals = Etavals[0:8] #index 8 is lagged middle cell itself, 9 is current value of middle cell
#                
#                I_xval = np.sum([np.cos(angles[i])*Etavals[i] for i in range(0,8)])
#                I_yval = np.sum([np.sin(angles[i])*Etavals[i] for i in range(0,8)])
#                
#                Eta_vectorstrength[la_ind,lo_ind,s]=np.sqrt(I_xval**2 + I_yval**2)
#                Eta_mainangle[la_ind,lo_ind,s] = np.arctan2(I_yval,I_xval)
#                           
#                #next retain the totol 4-d info from top 4 U+S sources
#                I_top4sources[la_ind,lo_ind,s]=info['I_top4']
#                
#                #next retain U1,U2,R,S, and Itotal from opposite direction sources
#                
#                #each has 4 different pieces...E-W, NE-SW, N-S, NW-SE
#                I_2directions = np.asarray(info['I_directions'])
#                U_sources = (np.asarray(info['UN'])+np.asarray(info['US']))/I_2directions
#                R_sources = np.asarray(info['R'])/I_2directions
#                S_sources = np.asarray(info['S'])/I_2directions
#                    
#                U_EW[la_ind,lo_ind,s]=U_sources[0]
#                U_NESW[la_ind,lo_ind,s]=U_sources[1]
#                U_NS[la_ind,lo_ind,s]=U_sources[2]
#                U_NWSE[la_ind,lo_ind,s]=U_sources[3]
#                
#                S_EW[la_ind,lo_ind,s]=S_sources[0]
#                S_NESW[la_ind,lo_ind,s]=S_sources[1]
#                S_NS[la_ind,lo_ind,s]=S_sources[2]
#                S_NWSE[la_ind,lo_ind,s]=S_sources[3]
#                
#                R_EW[la_ind,lo_ind,s]=R_sources[0]
#                R_NESW[la_ind,lo_ind,s]=R_sources[1]
#                R_NS[la_ind,lo_ind,s]=R_sources[2]
#                R_NWSE[la_ind,lo_ind,s]=R_sources[3]
#            
#                
#print('done with information theory analysis....')
#
#            
##%% loop over seasons, then info measure, then lat/lon to save maps
#                                        
#measures = [I_vectorstrength, I_mainangle, Eta_vectorstrength, Eta_mainangle,
#            I_top4sources, U_EW, U_NESW, U_NS, U_NWSE,S_EW,S_NESW,S_NS, S_NWSE,
#            R_EW, R_NESW, R_NS, R_NWSE] 
#measure_string =  ['Ivectorstrength', 'Imainangle', 'Etavectorstrength', 'Etamainangle',
#            'Itop4', 'U_EW', 'U_NESW', 'U_NS', 'U_NWSE','S_EW','S_NESW','S_NS', 'S_NWSE',
#            'R_EW', 'R_NESW', 'R_NS', 'R_NWSE']
#
#season_string = ['_DJF','_MAM','_JJA','_SON'] 
#
#filenames=[]
#allresults=[]
#
#for m_ind,m in enumerate(measures):
#    
#    print('now doing analysis for: '+ measure_string[m_ind])
#    
#    for s in range(0,4):
#        
#        data = m[:,:,s]
#        
#        filenames.append(measure_string[m_ind]+season_string[s]+'_avg')
#        allresults.append(data)
#                                
##%% now have a big list of filenames and allresults (list of lat-lon matrices)
#
##save in a pickle file
#
#site_filenamesave='results_directional_averages.pickle'                  
#f_myfile = open(site_filenamesave,'wb')
#pickle.dump(allresults, f_myfile)
#pickle.dump(filenames, f_myfile)
#pickle.dump(threshvals, f_myfile)
#f_myfile.close()   
#    
##%% save output as netcdf file...
#dsout = Dataset("SpatialInfo_Averages.nc4", "w", format="NETCDF4_CLASSIC")
#
#dsin = dataset
#
##Copy dimensions
#for dname, the_dim in dsin.dimensions.items():
#    print(dname, len(the_dim))
#    dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
#
#
## Copy variables
#for v_name, varin in dsin.variables.items():
#    outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
#    print(varin.datatype)
#    
#    # Copy variable attributes
#    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
#    
#    outVar[:] = varin[:]
#    
###add new variables, using same dimensions   
###vector angles (in radians)
#for a_ind,a in enumerate(allresults):    
#    X = dsout.createVariable(filenames[a_ind], np.float64, ('lat','lon',))  
#    X[:,:] = a        
#
#dsout.close()
    
    
    
    
    
    

#%% now want to look at "wind rose style plot" of I and Eta strength and direction


#%% saving output...netcdf files?
 
#input file
#dsin = dataset

#output file
#dsout = Dataset("SpatialInfo_5mmThreshold.nc4", "w", format="NETCDF4_CLASSIC")
#
##Copy dimensions
#for dname, the_dim in dsin.dimensions.items():
#    print(dname, len(the_dim))
#    dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
#
#
## Copy variables
#for v_name, varin in dsin.variables.items():
#    outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
#    print(varin.datatype)
#    
#    # Copy variable attributes
#    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
#    
#    outVar[:] = varin[:]
## close the output file
#    
#
##add new variables, using same dimensions
#    
##vector angles (in radians)
#Iangle_winter = dsout.createVariable('Iangle_winter', np.float64, ('lat','lon',))  
#Iangle_spring = dsout.createVariable('Iangle_spring', np.float64, ('lat','lon',)) 
#Iangle_summer = dsout.createVariable('Iangle_summer', np.float64, ('lat','lon',)) 
#Iangle_fall = dsout.createVariable('Iangle_fall', np.float64, ('lat','lon',)) 
#
#Ivectstrength_winter = dsout.createVariable('Ivectstrength_winter', np.float64, ('lat','lon',))  
#Ivectstrength_spring = dsout.createVariable('Ivectstrength_spring', np.float64, ('lat','lon',)) 
#Ivectstrength_summer = dsout.createVariable('Ivectstrength_summer', np.float64, ('lat','lon',)) 
#Ivectstrength_fall = dsout.createVariable('Ivectstrength_fall', np.float64, ('lat','lon',)) 
#
#Itop4_winter = dsout.createVariable('Itop4_winter', np.float64, ('lat','lon',))  
#Itop4_spring = dsout.createVariable('Itop4_spring', np.float64, ('lat','lon',)) 
#Itop4_summer = dsout.createVariable('Itop4_summer', np.float64, ('lat','lon',)) 
#Itop4_fall = dsout.createVariable('Itop4_fall', np.float64, ('lat','lon',))
#
#U_NS_winter = dsout.createVariable('U_NS_winter', np.float64, ('lat','lon',))  
#U_NS_spring = dsout.createVariable('U_NS_spring', np.float64, ('lat','lon',)) 
#U_NS_summer = dsout.createVariable('U_NS_summer', np.float64, ('lat','lon',)) 
#U_NS_fall = dsout.createVariable('U_NS_fall', np.float64, ('lat','lon',))
#
#U_EW_winter = dsout.createVariable('U_EW_winter', np.float64, ('lat','lon',))  
#U_EW_spring = dsout.createVariable('U_EW_spring', np.float64, ('lat','lon',)) 
#U_EW_summer = dsout.createVariable('U_EW_summer', np.float64, ('lat','lon',)) 
#U_EW_fall = dsout.createVariable('U_EW_fall', np.float64, ('lat','lon',))
#
#U_NESW_winter = dsout.createVariable('U_NESW_winter', np.float64, ('lat','lon',))  
#U_NESW_spring = dsout.createVariable('U_NESW_spring', np.float64, ('lat','lon',)) 
#U_NESW_summer = dsout.createVariable('U_NESW_summer', np.float64, ('lat','lon',)) 
#U_NESW_fall = dsout.createVariable('U_NESW_fall', np.float64, ('lat','lon',))
#
#U_NWSE_winter = dsout.createVariable('U_NWSE_winter', np.float64, ('lat','lon',))  
#U_NWSE_spring = dsout.createVariable('U_NWSE_spring', np.float64, ('lat','lon',)) 
#U_NWSE_summer = dsout.createVariable('U_NWSE_summer', np.float64, ('lat','lon',)) 
#U_NWSE_fall = dsout.createVariable('U_NWSE_fall', np.float64, ('lat','lon',))
#
#R_NS_winter = dsout.createVariable('R_NS_winter', np.float64, ('lat','lon',))  
#R_NS_spring = dsout.createVariable('R_NS_spring', np.float64, ('lat','lon',)) 
#R_NS_summer = dsout.createVariable('R_NS_summer', np.float64, ('lat','lon',)) 
#R_NS_fall = dsout.createVariable('R_NS_fall', np.float64, ('lat','lon',))
#
#R_EW_winter = dsout.createVariable('R_EW_winter', np.float64, ('lat','lon',))  
#R_EW_spring = dsout.createVariable('R_EW_spring', np.float64, ('lat','lon',)) 
#R_EW_summer = dsout.createVariable('R_EW_summer', np.float64, ('lat','lon',)) 
#R_EW_fall = dsout.createVariable('R_EW_fall', np.float64, ('lat','lon',))
#
#R_NESW_winter = dsout.createVariable('R_NESW_winter', np.float64, ('lat','lon',))  
#R_NESW_spring = dsout.createVariable('R_NESW_spring', np.float64, ('lat','lon',)) 
#R_NESW_summer = dsout.createVariable('R_NESW_summer', np.float64, ('lat','lon',)) 
#R_NESW_fall = dsout.createVariable('R_NESW_fall', np.float64, ('lat','lon',))
#
#R_NWSE_winter = dsout.createVariable('R_NWSE_winter', np.float64, ('lat','lon',))  
#R_NWSE_spring = dsout.createVariable('R_NWSE_spring', np.float64, ('lat','lon',)) 
#R_NWSE_summer = dsout.createVariable('R_NWSE_summer', np.float64, ('lat','lon',)) 
#R_NWSE_fall = dsout.createVariable('R_NWSE_fall', np.float64, ('lat','lon',))
#
#S_NS_winter = dsout.createVariable('S_NS_winter', np.float64, ('lat','lon',))  
#S_NS_spring = dsout.createVariable('S_NS_spring', np.float64, ('lat','lon',)) 
#S_NS_summer = dsout.createVariable('S_NS_summer', np.float64, ('lat','lon',)) 
#S_NS_fall = dsout.createVariable('S_NS_fall', np.float64, ('lat','lon',))
#
#S_EW_winter = dsout.createVariable('S_EW_winter', np.float64, ('lat','lon',))  
#S_EW_spring = dsout.createVariable('S_EW_spring', np.float64, ('lat','lon',)) 
#S_EW_summer = dsout.createVariable('S_EW_summer', np.float64, ('lat','lon',)) 
#S_EW_fall = dsout.createVariable('S_EW_fall', np.float64, ('lat','lon',))
#
#S_NESW_winter = dsout.createVariable('S_NESW_winter', np.float64, ('lat','lon',))  
#S_NESW_spring = dsout.createVariable('S_NESW_spring', np.float64, ('lat','lon',)) 
#S_NESW_summer = dsout.createVariable('S_NESW_summer', np.float64, ('lat','lon',)) 
#S_NESW_fall = dsout.createVariable('S_NESW_fall', np.float64, ('lat','lon',))
#
#S_NWSE_winter = dsout.createVariable('S_NWSE_winter', np.float64, ('lat','lon',))  
#S_NWSE_spring = dsout.createVariable('S_NWSE_spring', np.float64, ('lat','lon',)) 
#S_NWSE_summer = dsout.createVariable('S_NWSE_summer', np.float64, ('lat','lon',)) 
#S_NWSE_fall = dsout.createVariable('S_NWSE_fall', np.float64, ('lat','lon',))
#
#
#Iangle_winter[:,:] = I_mainangle[:,:,0]    
#Iangle_spring[:,:] = I_mainangle[:,:,1] 
#Iangle_summer[:,:] = I_mainangle[:,:,2] 
#Iangle_fall[:,:] = I_mainangle[:,:,3]
#
#Ivectstrength_winter[:,:] = I_vectorstrength[:,:,0]
#Ivectstrength_spring[:,:] =I_vectorstrength[:,:,1]
#Ivectstrength_summer[:,:] = I_vectorstrength[:,:,2]
#Ivectstrength_fall[:,:] =I_vectorstrength[:,:,3]
#
#Itop4_winter[:,:] = I_top4sources[:,:,0]
#Itop4_spring[:,:] = I_top4sources[:,:,1]
#Itop4_summer[:,:] =I_top4sources[:,:,2]
#Itop4_fall[:,:] = I_top4sources[:,:,3]
#
#U_NS_winter[:,:] = U_NS[:,:,0]
#U_NS_spring[:,:] = U_NS[:,:,1]
#U_NS_summer[:,:] = U_NS[:,:,2]
#U_NS_fall[:,:] = U_NS[:,:,3]
#
#U_EW_winter[:,:] =  U_EW[:,:,0]
#U_EW_spring[:,:] =U_EW[:,:,1]
#U_EW_summer[:,:] = U_EW[:,:,2]
#U_EW_fall[:,:] = U_EW[:,:,3]
#
#U_NESW_winter[:,:] = U_NESW[:,:,0]  
#U_NESW_spring[:,:] = U_NESW[:,:,1]
#U_NESW_summer[:,:] = U_NESW[:,:,2]
#U_NESW_fall[:,:] =U_NESW[:,:,3]
#
#U_NWSE_winter[:,:] = U_NWSE[:,:,0] 
#U_NWSE_spring[:,:] = U_NWSE[:,:,1]
#U_NWSE_summer[:,:] = U_NWSE[:,:,2]
#U_NWSE_fall[:,:] = U_NWSE[:,:,3]
#
#R_NS_winter[:,:] = R_NS[:,:,0]
#R_NS_spring[:,:] = R_NS[:,:,1]
#R_NS_summer[:,:] = R_NS[:,:,2]
#R_NS_fall[:,:] = R_NS[:,:,3]
#
#R_EW_winter[:,:] = R_EW[:,:,0] 
#R_EW_spring[:,:] =R_EW[:,:,1]
#R_EW_summer[:,:] = R_EW[:,:,2]
#R_EW_fall[:,:] =R_EW[:,:,3]
#
#R_NESW_winter[:,:] =R_NESW[:,:,0]
#R_NESW_spring[:,:] = R_NESW[:,:,1]
#R_NESW_summer[:,:] = R_NESW[:,:,2]
#R_NESW_fall[:,:] = R_NESW[:,:,3]
#
#R_NWSE_winter[:,:] = R_NWSE[:,:,0]
#R_NWSE_spring[:,:] =  R_NWSE[:,:,1]
#R_NWSE_summer[:,:] = R_NWSE[:,:,2]
#R_NWSE_fall[:,:] =R_NWSE[:,:,3]
#
#S_NS_winter[:,:] = S_NS[:,:,0]
#S_NS_spring[:,:] =  S_NS[:,:,1]
#S_NS_summer[:,:] = S_NS[:,:,2]
#S_NS_fall [:,:] =S_NS[:,:,3]
#
#S_EW_winter[:,:] = S_EW[:,:,0]
#S_EW_spring[:,:] = S_EW[:,:,1]
#S_EW_summer[:,:] = S_EW[:,:,2]
#S_EW_fall[:,:] =S_EW[:,:,3]
#
#S_NESW_winter[:,:] =S_NESW[:,:,0]
#S_NESW_spring[:,:] =S_NESW[:,:,1]
#S_NESW_summer[:,:] =S_NESW[:,:,2]
#S_NESW_fall[:,:] =S_NESW[:,:,3]
#
#S_NWSE_winter[:,:] =  S_NWSE[:,:,0]
#S_NWSE_spring[:,:] =S_NWSE[:,:,1]
#S_NWSE_summer[:,:] =S_NWSE[:,:,2]
#S_NWSE_fall[:,:] =S_NWSE[:,:,3]
#
#
#dsout.close()

#%%
   
            
#root_grp = Dataset('test5.nc', 'w', format='NETCDF4')
#root_grp.description = 'testing output of rainfall program'
#
## dimensions
#root_grp.createDimension('season', None)
#root_grp.createDimension('lon', len(lon))
#root_grp.createDimension('lat', len(lat))
#
## variables
#season = root_grp.createVariable('season', 'f8', ('season',))
#lon_vals = root_grp.createVariable('lon', 'f8', ('lon',))
#lat_vals = root_grp.createVariable('lat', 'f8', ('lat',))
#I_angle = root_grp.createVariable('I_angle', 'f8', ('lat','lon','season',))
#
#lon_vals=np.asarray(lon)
#lat_vals=np.asarray(lat)
#
## data
#for i in range(4):
#    season[i] = i
#    I_angle[:,:,i] = I_mainangle[:,:,i]
#
#root_grp.close()            
                
        
  
#%%      
#color map - make segmented colormap
#def make_segmented_cmap(): 
#    white = '#ffffff'
#    black = '#000000'
#    red = '#ff0000'
#    blue = '#0000ff'
#    anglemap = col.LinearSegmentedColormap.from_list(
#        'anglemap', [black, red, white, blue, black], N=256, gamma=1)
#    return anglemap        
#        
#N = 256
#cm = make_segmented_cmap()
#
#Iplot = I_mainangle[:,:,2]
#
#fig = plt.figure(1)
#fig.clf()
#plt.figure(figsize=(18,20))
#pic=plt.imshow(np.flipud(Iplot),cmap=cm)
#cbar = fig.colorbar(pic)
#cbar.ax.set_ylabel('Direction')
#
##%%
#
#Istrength = I_vectorstrength[:,:,0]
#
#fig = plt.figure(2)
#fig.clf()
#plt.figure(figsize=(12,15))
#plt.imshow(np.flipud(I_vectorstrength),vmin=0,vmax=0.4)




#%%



               
#now, separate these into watersheds - average watershed behavior                
                
                #maxinfoindex.append(np.argmax(vect))
                                             
#                localf = info['f']
#                locallogf = info['logf']
#                localplogf = info['plogf']
#                print(np.sum(localplogf))
#                
#                LocalVals[la_ind, lo_ind, p_ind, :,:,:,:,:,:]=localplogf
#                Localf[la_ind, lo_ind, p_ind, :,:,:,:,:,:]=localf
#                LocalFracVals[la_ind, lo_ind, p_ind, :,:,:,:,:,:]=localplogf/info['Itotal']
#                
#                
#                
#                posvals = np.zeros(np.shape(localplogf))
#                negvals = np.zeros(np.shape(localplogf))
#                
#                vals = np.reshape(localplogf,(1,64))
#                PosTot=np.sum(vals[vals>0])
#                NegTot=np.sum(vals[vals<0])
#                
#  
#                TotalPosVals[la_ind,lo_ind,p_ind]=PosTot
#                TotalNegVals[la_ind,lo_ind,p_ind]=NegTot
        
 



               
#%%
#n=0
#row = 0
#
#fig1 = plt.figure(figsize=(12,12))
#
#
#for i in range(0,2):
#    for j in range(0,2):
#        for k in range(0,2):
#            for m in range(0,2):
#                fig1.add_subplot(4,4,n+1)
#                test = LocalFracVals[:,:,0,i,j,k,m]                
#                plt.imshow(np.flipud(test),vmin=-2, vmax=2,cmap='bwr')
#    
#                plt.colorbar()
#                plt.title('I local (' + str(i)+str(j)+ str(k)+ str(m)+ ') /Itotal')
#                
#                n+=1
#         
# 
##%%           
#            
#fig2 = plt.figure(figsize=(10,10))            
#n =0           
#for i in range(0,2):
#    for j in range(0,2):
#        for k in range(0,2):            
#            for m in range(0,2):
#                fig2.add_subplot(4,4,n+1)
#                test = Localf[:,:,0,i,j,k,m]                
#                plt.imshow(np.flipud(test),vmin=0, vmax=2,cmap='bwr')
#    
#                plt.colorbar()
#                plt.title('p('+ str(i)+str(j)+ str(k)+ str(m)+ ') / p('+str(i)+str(j)+str(k)+')p('+str(m)+')')
#                
#                n+=1

            
  #%%  
#  
#fig3 = plt.figure(figsize=(15,30))            
#nn =0           
#for i in range(0,2):
#    for j in range(0,2):
#        for k in range(0,2):            
#            for m in range(0,2):
#                for n in range(0,2):
#                    for o in range(0,2):
#                        fig3.add_subplot(16,4,nn+1)
#                        test = LocalVals[:,:,0,i,j,k,m,n,o] 
#                    
#                           
#                        plt.imshow(np.flipud(test),vmin=-.2, vmax=.2,cmap='bwr')
#            
#                        #plt.colorbar()
#                        plt.title('p('+ str(i)+str(j)+ str(k)+ str(m)+str(n)+str(o)+')log[p('+ str(i)+str(j)+ str(k)+str(m)+str(n) +
#                                    str(o)+ ') / p('+str(i)+str(j)+str(k)+str(m)+str(n)+')p('+str(o)+')]')
#            
#                        nn+=1
#
#    
##%%  
#
#
#fig4 = plt.figure(figsize=(25,15))
#fig4.add_subplot(4,1,1)
#plt.imshow(np.flipud(TotalPosVals[:,:,0]),vmin=0,vmax=.6)
#plt.colorbar()
#plt.title('Positive Part of Total I')
#fig4.add_subplot(4,1,2)
#plt.imshow(np.flipud(-TotalNegVals[:,:,0]),vmin=0,vmax=.6)
#plt.colorbar()
#plt.title('Negative Part of Total I')
#fig4.add_subplot(4,1,3)
#ratio = TotalPosVals[:,:,0]/-TotalNegVals[:,:,0]
#plt.imshow(np.flipud(ratio),vmin=1,vmax=3)
#plt.colorbar()
#plt.title('Ratio of Positive/Negative Parts')
#totalI = TotalPosVals[:,:,0]+TotalNegVals[:,:,0]        
#fig4.add_subplot(4,1,4)
#plt.imshow(np.flipud(totalI),vmin=0,vmax=.4)
#plt.colorbar() 
#plt.title('Total I')        