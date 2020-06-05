#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 10:31:22 2018

This program computes multidimensional lagged MI between lagged and current rainfall
histories -

Rainfall from CPC gage-based gridded data, converted to binary according to 
%-ile threshold (after values below 0.3 have been removed)

consider lagged mutual information between neighboring grid cells (8 directions)
resolve vector of information values to obtain "dominant direction"

%dominant direction also associated with dominant strength
%but also have strength of information from all 8 directions, and self (or choose top 3-5 directions?)

%Addition in Feb 2019: import file of climate indices for seasons - 
%look for correlations between climate indices and information values
%for 3 year windows

%Feb 2019: split into two codes - one for overall, one for trends and correlations
%Feb 2019: split into another segment - one code just to determine lower threshold values

%THIS code detects thresholds and saves as a pickle file
% threhsolds are based on the maximally infrormative rainfall magnitude

@author: allisongoodwell
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import scipy.stats as stats
import pickle

import time


np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

N = 2 #number of bins (2 for binary pdf)
years = range(1948,2018)
nyrs = len(years)

#function to compute mutual information AND statistical significance, 2d case only
def compute_info_measures(Tuple):
    
    ntests=10
    sig_level = 3

    pdf,edges = np.histogramdd(Tuple,bins=2)
    pdf = pdf/np.sum(pdf)
    
    source_dim = 0;
    tar_dim = 1;
    
    m_sources = np.sum(pdf,axis=tar_dim) #should be dim-1 dimension pdf
    m_target = np.sum(pdf,axis=source_dim) #should be 1D pdf
    
    Hx = -m_target[0]*np.log2(m_target[0]) - m_target[1]*np.log2(m_target[1])
   
    #pad m_sources with 1 for last index
    m_sources = np.expand_dims(m_sources, axis=np.size(np.shape(m_sources)))
         
    p_independent = m_sources * m_target
    f_ijk = pdf/p_independent
    
    logf = np.where(p_independent != 0, np.log2(f_ijk), 0)  
    plogf = pdf*logf

    Ival = np.nansum(plogf)  #return normalized value of I
    
    #shuffle all dimensions except last one (the target)
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
        f_ijk = pdf/p_independent
        
        logf = np.where(p_independent != 0, np.log2(f_ijk), 0)  
        plogf = pdf*logf
        
        Itot_vector[j] = np.nansum(plogf)
    
    I_mean=np.average(Itot_vector)
    I_std=np.std(Itot_vector)  
    Icrit = I_mean + sig_level*I_std
    
    if Icrit>Ival:
        Ival=0
    else:
        Ival = Ival/Hx
        
    return Ival

#%%

ppt_list=[]
for y in years:
    
    string = 'C:\\Users\\goodwela\\Dropbox\\UCDenver\\rainfall_research\\CPC_raingage_gridded_dailydata\\original_precip_data\\precip.V1.0.'+str(y)+'.nc'
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

#based on maximum Lagged information
thresh_vect =[0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9]

threshvals = np.zeros((np.size(lat),np.size(lon),4)) 
avg_dailymagnitude = np.zeros((np.size(lat),np.size(lon),4)) 
avg_rainymagnitude = np.zeros((np.size(lat),np.size(lon),4)) 

for la_ind,la in enumerate(lat):
    
    if np.mod(la_ind,10)==0:
        print(la)
    
    if la_ind<1:
        continue
    if la == np.max(lat):
        continue
                
    for lo_ind,lo in enumerate(lon):
                
        if lo_ind<1:
            continue
        if lo == np.max(lon):
            continue
            
  
        #for each latitude and longitude, extract precip data  
        for s in range(0,4):    
                          
            ppt_test = ppt_list[0][0,la_ind,lo_ind]

            if np.isnan(ppt_test)==True:
           
                continue
            
            ppt_list_mini = [ppt_list[y_ind][days[s],la_ind,lo_ind] for y_ind,y in enumerate(years)]
            PPTvect = np.concatenate(ppt_list_mini)
            
            
            #indices of season breaks (new season stars, old season ends)
            cutoff_starts = [y_ind*len(days[s])-1 for y_ind,y in enumerate(years)]
                   
            if np.isnan(PPTvect[0])==True:
                #print('skipping!')
                continue
            #PPTvect is all precip for a given coordinate            
            PPTvect[np.isnan(PPTvect)]=0
            
            Xtar = PPTvect[1:]
            Xlag = PPTvect[:-1]
            #remove season breaks from Xtar and Xlag (need to find indices)
            Xtar = [Xtar[i] for i in range(0,len(Xtar)) if i not in cutoff_starts]
            Xlag = [Xlag[i] for i in range(0,len(Xlag)) if i not in cutoff_starts]
            

            #go through threshold values, find lagged MI
            MI_thresh=[]
            Tuple_list=[]
            
            for t in thresh_vect:

                #create binary variables for the threshold value t
                Xtar_binary = np.asarray(Xtar)      
                Xlag_binary = np.asarray(Xlag)
                Xtar_binary[Xtar_binary<t]=0
                Xtar_binary[Xtar_binary>=t]=1
                Xlag_binary[Xlag_binary<t]=0
                Xlag_binary[Xlag_binary>=t]=1
                          
                Tuple = np.vstack((Xlag_binary,Xtar_binary))            
                Tuple = np.transpose(Tuple)
                
                Tuple_list.append(Tuple)
                  
            t1 = time.clock()
            
            print('going into I calcs with '+ str(len(Tuple_list))+ ' for latitude '+ str(la) +'and long ' +str(lo))  
            Ivals=[]
            for T in Tuple_list:
                I = compute_info_measures(T)
                Ivals.append(I)
            
            
            t2 = time.clock() 
            print('time elapsed= '+ str(t2-t1))
            
            #print(Ivals)
                        
            threshvals[la_ind,lo_ind,s]=thresh_vect[np.argmax(Ivals)]
            
            newPPTvect = PPTvect[PPTvect>thresh_vect[np.argmax(Ivals)]]
            avg_dailymagnitude[la_ind,lo_ind,s]=np.mean(PPTvect)
            avg_rainymagnitude[la_ind,lo_ind,s]=np.mean(newPPTvect)
                
print('threshold values determined....')

#%%
#save thresholdsa as a netcdf file and a pickle file
dsout = Dataset("Spatial_Thresholds_Avgs_mm.nc4", "w", format="NETCDF4_CLASSIC")

dsin = dataset

#Copy dimensions
for dname, the_dim in dsin.dimensions.items():
    print(dname, len(the_dim))
    dsout.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)


# Copy variables
for v_name, varin in dsin.variables.items():
    outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
    print(varin.datatype)
    
    # Copy variable attributes
    outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    
    outVar[:] = varin[:]
    
##add new variables, using same dimensions   
##vector angles (in radians)
filenames = ['Thresholds_DJF','Thresholds_MAM','Thresholds_JJA','Thresholds_SON',
             'AvgMag_DJF','AvgMag_MAM','AvgMag_JJA','AvgMag_SON',
             'AvgRainyMag_DJF','AvgRainyMag_MAM','AvgRainyMag_JJA','AvgRainyMag_SON',]
for s in range(0,4):  
    X = dsout.createVariable(filenames[s], np.float64, ('lat','lon',))  
    X[:,:] = threshvals[:,:,s]        
for s in range(0,4):  
    X = dsout.createVariable(filenames[s+4], np.float64, ('lat','lon',))  
    X[:,:] = avg_dailymagnitude[:,:,s]     
for s in range(0,4):  
    X = dsout.createVariable(filenames[s+8], np.float64, ('lat','lon',))  
    X[:,:] = avg_rainymagnitude[:,:,s]     

dsout.close()

site_filenamesave='Spatial_thresholds_mm.pickle'                  
f_myfile = open(site_filenamesave,'wb')
pickle.dump(threshvals, f_myfile)
f_myfile.close() 

print('threshold values saved....')
    
    
    #%%
    #Iplot = threshvals[:,:,3]
    #fig = plt.figure(1)
    #plt.figure(figsize=(12,15))
    #pic=plt.imshow(np.flipud(Iplot))        
    
