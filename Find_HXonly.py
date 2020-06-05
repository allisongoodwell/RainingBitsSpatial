#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 10:31:22 2018

This program computes ONLY the entropy for binary precip, entire time period

%This one takes threshold values as inputs

@author: allisongoodwell
"""

from netCDF4 import Dataset
import numpy as np
import pickle


np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})


lags = [1]
N = 2 #number of bins (2 for binary pdf)
years = range(1948,2018)
nyrs = len(years)

#function to compute seasonal HX
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
        
        p_rain = np.sum(Xtar)/len(Xtar)
        
        Hx = -p_rain*np.log2(p_rain) - (1-p_rain)*np.log2(1-p_rain)
                                      
        infodict = {'H_x':Hx}
        
        
        info_all.append(infodict) #list of dictionaries (one for each season)
     
    return info_all

#%%

ppt_list=[]
for y in years:
    
    string = 'C:\\Users\\goodwela\\Dropbox\\UCDenver\\research\\rainfall_research\\CPC_raingage_gridded_dailydata\\original_precip_data\\precip.V1.0.'+str(y)+'.nc'
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

Hxvals = np.empty((np.size(lat),np.size(lon),4))*np.nan

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
               
        ppt_test = ppt_list[0][0,la_ind,lo_ind]
        
        if lo_ind<1 or lo==np.max(lon) or np.isnan(ppt_test)==True:
            continue
        
        print('longitude = '+str(lo) + ' latitude = ' +str(la))
        
        
        ppt_data = [ppt_list[y_ind][:,la_ind,lo_ind] for y_ind,y in enumerate(years)]
               
    
        for s in range(0,4):
            
            #for each latitude and longitude pair, find 8 nearest neighbors                 
 
            PPTvect = np.concatenate(ppt_data)
            
            cutoff_starts = [y*len(days[s])-1 for y in range(1,len(years))]
            
            if np.isnan(PPTvect[0])==True:
                #print('skipping!')
                continue
            #PPTvect is all precip for a given coordinate, and 8 surrounding coordinates
            # size time length x latitudes x longitudes
            
            PPTvect[np.isnan(PPTvect)]=0
            
            #get binary vector of all rainfall values
            pptvect_binary = np.zeros(np.shape(PPTvect))
            
            pptvect_binary[PPTvect>threshvals[la_ind,lo_ind,s]]=1         
        
            Xtar = pptvect_binary       
            Xtar = np.delete(Xtar,cutoff_starts,axis=0)
            
            p_rain = np.sum(Xtar)/len(Xtar)
            
            Hx = -p_rain*np.log2(p_rain) - (1-p_rain)*np.log2(1-p_rain)
                                          
            Hxvals[la_ind,lo_ind,s]=Hx
            
            
site_filenamesave = 'overall_Hx_values.pickle'
f_myfile = open(site_filenamesave,'wb')
pickle.dump(Hxvals,f_myfile)
f_myfile.close()   

print('saving results pickle file')
        
    