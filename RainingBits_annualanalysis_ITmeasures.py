#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated May 2020

This program computes multidimensional lagged MI between lagged and current rainfall
histories, on an annual window, for four seasons
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
pickle files for each row of latitude, in folder called TrendResults

Notes:
*this code runs compute_info_measures function in parallel for each row of latitude
must choose number of processors in pool
*next run code CPC_trendscorrelations.py to analyze trends over time and correlations with climate indices

@author: allisongoodwell
"""

from netCDF4 import Dataset
import numpy as np
import pickle
from multiprocessing import Pool
import time
import gc


precip_folder='C:\\Users\\goodwela\\Dropbox\\UCDenver\\research\\rainfall_research\\CPC_raingage_gridded_dailydata\\original_precip_data'
nProcessors = 14


np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

maxdim = 9 #9 lagged neighbors (including self) and 1 current target at center
N = 2 #number of bins (2 for binary pdf)
years = range(1948,2018)


#function to compute statistical significance, shuffled surrogates
def compute_stat_sig(Tuple):
    
    #shuffle all dimensions except last one (the target)
    ntests=100
    sig_level = 3

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
    threshvals = input_list[2]
    days = input_list[3]
    
    info_all=[]
    
    y_step=1 #aggregate years (e.g. 1 year vs 3 year averages vs 5 year)
    years_int = range(1950,2018,y_step)
        
    for y_ind,y in enumerate(years_int):
    
        y_start = y_ind*y_step
        year_vals = range(y_start,(y_start+y_step))
        
        info_years=[]       
        for s in range(0,4):
            
            #for each latitude and longitude pair, find 8 nearest neighbors                 
            ppt_list_mini = [ppt_list[x][days[s],:,:] for x in year_vals]
            PPTvect = np.concatenate(ppt_list_mini)
            
            cutoff_starts = [y*len(days[s])-1 for y in range(1,y_step)]
            
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
            
            #find critical I based on self-dependency for center (and use these for neighbors)
            Icrit=0;
            if np.sum(Tuple[:,8])>=1:
                smallTuple = Tuple[:,[8,-1]]           
                Icrit = compute_stat_sig(smallTuple)

            
            
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
                
                I_here = I_2D[-1] #lagged info between center cell and itself is last index of I_2D
                #check for statistical significance of Ival (but only if I_here is 0)
                if I_here ==0:
                    Icrit = compute_stat_sig(Tuple3d)
                    if Icrit > Ival:
                        Ival=0
                        Eta.append(0)
                    else:
                        Eta.append(np.divide(Ival-I_here,Ival))
                
                I_3D.append(Ival)
                
                #now break this down into U_here, U_there, R and S components
                #determine eta = U_there + S / Itot as the information proportion above that from center itself
                
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
                        Itotal=0
                    
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
            info_years.append(infodict) #list of dictionaries (one for each season)
        
        info_all.append(info_years) #(append for each year)
     
    return info_all

#%% main code - runs function above in parallel for grid cells in each row of longitude
    
if __name__ == '__main__':
    ppt_list=[]
    for y in years:
        
        string = precip_folder + '\\precip.V1.0.'+str(y)+ '.nc'
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
    
    thresh_file='Spatial_thresholds_mm.pickle'                  
    f_myfile = open(thresh_file,'rb')
    threshvals = pickle.load(f_myfile)
    f_myfile.close() 
           
    #%% 
    
    maxinfoindex=[]
    
    maxinfoindex=[]
    AllResults=[]
    Long_Lists=[]
    Lat=[]
        
    for la_ind,la in enumerate(lat):
             
        if la_ind<1 or la == np.max(lat):
            continue
        
        
        pdf_long_list=[]
        ppt_long_list=[]
        
    
        data_list=[]
        lon_list=[]
        lat_list=[] 
        item_list=[]    
        for lo_ind,lo in enumerate(lon):
            
            #t1 = time.clock()

            
            if lo_ind<1 or lo==np.max(lon):
                continue
            
            #goal to make this parallel -- save ppt_list as set of lists for each lat, lon pair
            ppt_test = ppt_list[0][0,la_ind,lo_ind]
            #print(ppt_test)
            if np.isnan(ppt_test)==True:
                #t2 = time.clock()
                #print('skipping! time elapsed'+str(t2-t1))                
                continue
             
                
            print('longitude = '+str(lo) + ' latitude = ' +str(la))
            
            ppt_data = [ppt_list[y_ind][:,(la_ind-1):(la_ind+2),(lo_ind-1):(lo_ind+2)] for y_ind,y in enumerate(years)]
                       
            item_list = [ppt_data, years, threshvals[la_ind,lo_ind,:],days]
                   
            data_list.append(item_list)
            lon_list.append(lo)
            #function that takes an element from data_list

        gc.collect()    
        t1 = time.clock()
        print('going into metrics with '+ str(len(lon_list))+ ' for latitude '+ str(la))    
        p = Pool(nProcessors)
        results = p.map(compute_info_measures,data_list) #pdfs should be a list of lists
        p.close
        p.join
      
        t2 = time.clock()  
        print('time elapsed for latitude row = '+ str(t2-t1))    
        
        gc.collect()
    
        if len(lon_list)>0:    
            site_filenamesave = 'TrendResults\\spatial_results_trends_'+str(la)+'.pickle'
            f_myfile = open(site_filenamesave,'wb')
            pickle.dump(results,f_myfile)
            pickle.dump(lon_list,f_myfile)
            pickle.dump(la,f_myfile)
            f_myfile.close()
             
        
#        AllResults.append(results)   
#        Long_Lists.append(lon_list)
#        Lat.append(la)
#        
#        
#    site_filenamesave = 'spatial_results_trends.pickle'
#    f_myfile = open(site_filenamesave,'wb')
#    pickle.dump(AllResults,f_myfile)
#    pickle.dump(Long_Lists,f_myfile)
#    pickle.dump(Lat,f_myfile)
#    f_myfile.close()   
#    
#    print('saving results pickle file')
    
    

            
