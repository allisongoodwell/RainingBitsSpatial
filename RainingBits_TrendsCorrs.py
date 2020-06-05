#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 10:31:22 2018

This program uses results from CPC_annual_ITanalysis_051019.py
(multiple pickle files located in TrendResults folder)

takes annual IT measures and detects trends in several selected metrics

inputs:
    one CPC precipitation file (e.g. precip.V1.0.1948.nc here)
    folder of Trend Results 'TrendResults/spatial_results_trends_'+str(la)+'.pickle'
    where str(la) is a latitude value (1/8 degree resolution)
    climateinds_10.pickle, which contains seasonally averaged climate indices
    
    
outputs:
    file of all results: SpatialInfoTrendResults.pickle                  
    netcdf file of trends for different measures: SpatialInfoTrendMaps.nc4
    excel file of averages for each metric: SpatialInfoTrendAverages.xlsx



@author: allisongoodwell
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import chi2
import pickle

#from mpl_toolkits.basemap import Basemap

np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

y_step=1 #aggregate years (e.g. 1 year vs 3 year averages vs 5 year)
 
years = range(1950,2018,y_step)
nyrs = len(years)
    
precip_filename = 'precip.V1.0.1948.nc' #only need for lat and lon values
dataset = Dataset(precip_filename)
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]


#vector resolution: find x and y components of directions
angles = np.asfarray([90, 135, 180, 225, 270, 315, 0, 45])
angles = angles*np.pi/180



#%%

#here, have many pickle files to load in with results, then need to use to find trends

I_2dvals = np.empty((np.size(lat),np.size(lon),nyrs,4,9))*np.nan
I_3dvals = np.empty((np.size(lat),np.size(lon),nyrs,4,9))*np.nan
Eta_vals = np.empty((np.size(lat),np.size(lon),nyrs,4,9))*np.nan
I_mainangle = np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
I_vectorstrength = np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
Eta_mainangle = np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
Eta_vectorstrength = np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan

I_top4sources = np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan

U_EW= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
U_NS= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
U_NESW= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
U_NWSE= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan

S_EW= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
S_NS= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
S_NESW= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
S_NWSE= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan

R_EW= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
R_NS= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
R_NESW= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan
R_NWSE= np.empty((np.size(lat),np.size(lon),nyrs,4))*np.nan


for la_ind,la in enumerate(lat):
    
    if np.mod(la_ind,10)==0:
        print('latitude = '+str(la))
    
    if la<25.2 or la>49.3:
        continue
      
    #site_filenamesave = 'TrendResults\\spatial_results_trends_'+str(la)+'.pickle'
    site_filenamesave = 'TrendResults/spatial_results_trends_'+str(la)+'.pickle'
    f_myfile = open(site_filenamesave,'rb')
    results = pickle.load(f_myfile)
    lon_list = pickle.load(f_myfile)
    latval = pickle.load(f_myfile)
    f_myfile.close()
    
    #now want to go through longitudes and fill in matrices with data
    i=0
    for lo_ind, lo in enumerate(lon):
        if i<len(lon_list) and lo == lon_list[i]:  #if longitude matches that of the next list in the row of results
            
            res = results[i]
            for y in range(nyrs):
                for s in range(4):

                    if s >= len(res[y]):
                        continue
                    
                    info = res[y][s]
                    
                    Ivals = info['I_2d']
                    Ivals = Ivals[0:8] #index 8 is lagged middle cell itself, 9 is current value of middle cell
                    
                    I_yval = np.sum([np.cos(angles[i])*Ivals[i] for i in range(0,8)])
                    I_xval = np.sum([np.sin(angles[i])*Ivals[i] for i in range(0,8)])
                    
                    I_vectorstrength[la_ind,lo_ind,y,s]=np.sqrt(I_xval**2 + I_yval**2)
                    if I_xval>0:
                        I_mainangle[la_ind,lo_ind,y,s] = np.arctan2(I_xval,I_yval)
                    else:
                        I_mainangle[la_ind,lo_ind,y,s] = 2*np.pi + np.arctan2(I_xval,I_yval)
                    
                    I_mainangle[la_ind,lo_ind,y,s]=I_mainangle[la_ind,lo_ind,y,s]*180/np.pi
                    
                    #print('Ivect strength = '+str(I_vectorstrength[la_ind,lo_ind,s]))
                    #compute primary direction (angle) and strength of info to center based on TE metrics
                    Etavals = info['Eta']
                    Etavals = Etavals[0:8] #index 8 is lagged middle cell itself, 9 is current value of middle cell
                    
                    I_xval = np.sum([np.cos(angles[i])*Etavals[i] for i in range(0,8)])
                    I_yval = np.sum([np.sin(angles[i])*Etavals[i] for i in range(0,8)])
                    
                    Eta_vectorstrength[la_ind,lo_ind,y,s]=np.sqrt(I_xval**2 + I_yval**2)
                    if I_xval>0:
                        Eta_mainangle[la_ind,lo_ind,y,s] = np.arctan2(I_xval,I_yval)
                    else:
                        Eta_mainangle[la_ind,lo_ind,y,s] = 2*np.pi + np.arctan2(I_xval,I_yval)
                    
                    Eta_mainangle[la_ind,lo_ind,y,s]=Eta_mainangle[la_ind,lo_ind,y,s]*180/np.pi
                               
                    #next retain the totol 4-d info from top 4 U+S sources
                    I_top4sources[la_ind,lo_ind,y,s]=info['I_top4']
                    
                    #next retain U1,U2,R,S, and Itotal from opposite direction sources
                    
                    #each has 4 different pieces...E-W, NE-SW, N-S, NW-SE
                    I_2directions = np.asarray(info['I_directions'])
                    U_sources = (np.asarray(info['UN'])+np.asarray(info['US']))/I_2directions
                    R_sources = np.asarray(info['R'])/I_2directions
                    S_sources = np.asarray(info['S'])/I_2directions
                        
                    U_EW[la_ind,lo_ind,y,s]=U_sources[0]        
                    U_NS[la_ind,lo_ind,y,s]=U_sources[2]
                    
                    S_EW[la_ind,lo_ind,y,s]=S_sources[0]
                    S_NS[la_ind,lo_ind,y,s]=S_sources[2]
                   
                    R_EW[la_ind,lo_ind,y,s]=R_sources[0]
                    R_NS[la_ind,lo_ind,y,s]=R_sources[2]
                    
            i=i+1 #increment i to look for next lon value
          
#%% now, want to find correlation between climate indices and information values....
                
#need to load pickle file
string = 'climateinds_10.pickle'
fileObject = open(string,'rb')
b=pickle.load(fileObject)


Indices_3yr=np.zeros((4,nyrs,4)) #index, year, season
Corr = np.zeros((np.size(lat),np.size(lon),4,4)) #lat,lon,index,season

for i_ind,i in enumerate([b['NAO'],b['EPNP'],b['AMO'],b['PDO']]): 

    for y in range(nyrs):
        y_start = y*y_step
        year_vals = range(y_start,(y_start+y_step)) 
               
        for s in range(0,4):        
            Indices_3yr[i_ind, y,s]=np.mean(i[year_vals,s])
  
          

#%%
#loop over seasons, then info measure, then climate index, then lat/lon
                                   
measures = [I_vectorstrength, I_mainangle, Eta_vectorstrength, Eta_mainangle,
            I_top4sources, R_EW] 
measure_string =  ['Ivectorstrength', 'Imainangle', 'Etavectorstrength', 'Etamainangle',
            'Itop4', 'R_EW']
season_string = ['_DJF','_MAM','_JJA','_SON'] 

aspect_string = ['_trend','_corrNAO','_corrEPNP','_corrAMO','_corrPDO'] 

filenames=[]
allresults=[]
allmeanvals=[]

for m_ind,m in enumerate(measures):
    
    
    print('now doing analysis for: '+ measure_string[m_ind])
    
    for s in range(0,4):
        print(s)
        data = m[:,:,:,s]
        
        #find trends over time (must go with lat lon loops)
        Trend = np.zeros((np.size(lat),np.size(lon)))
        MeanVal = np.zeros((np.size(lat),np.size(lon)))
        n_slopes=0
        for la_ind,la in enumerate(lat):
                
                                
                if la_ind<1:
                    continue
                if la == np.max(lat):
                    continue
                        
                for lo_ind,lo in enumerate(lon):
                    
                    InfoVect = m[la_ind,lo_ind,:,s]
                    
                    #find average value of InfoVect
                    mean_infovect = np.nanmean(InfoVect)
                    
                    if np.isnan(InfoVect[0])==True:
                        continue
                 
                    
                    if m_ind==1 or m_ind==3: 
                        #for circular data - use modified linear regression 
                        RadianData = InfoVect * (np.pi/180)
                        R_x = np.mean(np.cos(RadianData))
                        R_y = np.mean(np.sin(RadianData))
                        
                    
                        if R_x>0:
                            mainangle = np.arctan2(R_x,R_y)
                        else:
                            mainangle = 2*np.pi + np.arctan2(R_x,R_y)
                    
                        mainangle=mainangle*180/np.pi
                        
                        print('main angle:')
                        print(mainangle)
                        
                        #redefine mean as the circular mean...
                        mean_infovect = mainangle
                        
                        #want to shift the values so that R_mean is 180 (or pi)
                        InfoVect = InfoVect-mainangle
                    
                    
                    #print(str(la) + ' ' +str(lo))
                    
#                    SenSlope = stats.theilslopes(InfoVect,years, 0.95)        
#                    sign_slope = np.sign(SenSlope[2])*np.sign(SenSlope[3])
                    
                    lsq, intercept, rval, pval,err = stats.linregress(years,InfoVect)
                    if pval<0.05:
                        Trend[la_ind,lo_ind]=lsq
                        MeanVal[la_ind,lo_ind]=mean_infovect
                    
#                    if sign_slope<0:
#                        SenSlope=0
#                    else:
#                        SenSlope=SenSlope[0]
#                        n_slopes=n_slopes+1
#                        #print(n_slopes)
#                    Trend[la_ind,lo_ind]=SenSlope

 
        #want to save Trend
        filenames.append(measure_string[m_ind]+aspect_string[0]+season_string[s])
        allresults.append(Trend)
        allmeanvals.append(MeanVal)
        
        #third: find correlations with climate indices

        for a_ind,a in enumerate([b['NAO'],b['EPNP'],b['AMO'],b['PDO']]):  
            #also must go with lat, lon loops to find correlations...
             #now need to go through every point and find the correlation
            IndexVect = Indices_3yr[a_ind,:,s]
            Corr = np.zeros((np.size(lat),np.size(lon)))
            
            for la_ind,la in enumerate(lat):
                                
                if la_ind<1:
                    continue
                if la == np.max(lat):
                    continue
                        
                for lo_ind,lo in enumerate(lon):
                    
                    InfoVect = m[la_ind,lo_ind,:,s]
                    
                    if np.isnan(InfoVect[0])==True:
                        continue
                    
                    #find circular correlation for directional variables
                    if m_ind==1 or m_ind==3:
                        RadianData = InfoVect * (np.pi/180)
                        R_indexinfox = np.corrcoef(IndexVect,np.cos(RadianData))[0,1]
                        R_indexinfoy = np.corrcoef(IndexVect,np.sin(RadianData))[0,1]
                        R_xy = np.corrcoef(np.cos(RadianData),np.sin(RadianData))[0,1]
                        
                        R = np.sqrt(((R_indexinfox**2)+(R_indexinfoy**2)-2*R_indexinfox*R_indexinfoy*R_xy)/(1-R_xy**2))

                    else: #just find linear correlation coefficience
                        R = np.corrcoef(IndexVect,InfoVect)[0,1]
                        
                    teststat = np.abs(R**2*len(IndexVect))
                    pval = chi2.sf(teststat,2) 
#                    if m_ind==1:
#                        print(R_indexinfox, R_indexinfoy, R_xy)
#                        print(R,pval)

                    if pval < 0.05:
                        
                        Corr[la_ind,lo_ind] = R
             
            #want to save Corr    
            filenames.append(measure_string[m_ind]+aspect_string[a_ind+1]+season_string[s])
            allresults.append(Corr)
            allmeanvals.append(MeanVal)
            
            
#%% now have a big list of filenames and allresults (list of lat-lon matrices)

#need to figure out how to save this...
#save in a pickle file
import pickle
#import imageio

site_filenamesave='SpatialInfoTrendResults.pickle'                  
f_myfile = open(site_filenamesave,'wb')
pickle.dump(allresults, f_myfile)
pickle.dump(allmeanvals, f_myfile)
pickle.dump(filenames, f_myfile)
f_myfile.close()   

#imageio.imwrite(filenames[0]+'.tif',allresults[0])        
#%%
#save output as netcdf file...
dsout = Dataset("SpatialInfoTrendMaps.nc4", "w", format="NETCDF4_CLASSIC")

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
# close the output file 
    
##add new variables, using same dimensions
#    
##vector angles (in radians)
for a_ind,a in enumerate(allresults):    
    X = dsout.createVariable(filenames[a_ind], np.float64, ('lat','lon',))  
    X[:,:] = a        

dsout.close()

#%% want to find averages for each thing...
#for each item in all_results: find number of sig values, and average value
numpos=[]
numneg=[]
avgpos=[]
avgneg=[]
avgtrend=[]
avgvaluepostrend=[]
avgvaluenegtrend=[]
avgpospercent=[]
avgnegpercent=[]

for a_ind,a in enumerate(allresults):
    print(filenames[a_ind])
    
    values = allmeanvals[a_ind]
    
    numpos.append(len([x for x in np.reshape(a,(np.size(a),)) if x>0]))
    numneg.append(len([x for x in np.reshape(a,(np.size(a),)) if x<0]))

    avgpos.append(np.average([x for x in np.reshape(a,(np.size(a),)) if x>0]))
    avgneg.append(np.average([x for x in np.reshape(a,(np.size(a),)) if x<0]))
    
    avgtrend.append(np.average([x for x in np.reshape(a,(np.size(a),)) if ~np.isnan(x)]))
    
    avgvaluepostrend.append(np.average([x for x in np.reshape(values,(np.size(values),)) if x>0]))
    avgvaluenegtrend.append(np.average([x for x in np.reshape(values,(np.size(values),)) if x<0]))
    
    percents = 100*(a*32)/np.abs(values)
    
    avgpospercent.append(np.average([x for x in np.reshape(percents,(np.size(percents),)) if x>0]))
    avgnegpercent.append(np.average([x for x in np.reshape(percents,(np.size(percents),)) if x<0]))

import pandas as pd
df = pd.DataFrame({'name':filenames,'avgtrend':avgtrend,'avgposvalue':avgvaluepostrend,'numpos':numpos,'avgpos':avgpos,
                   'pospercent':avgpospercent,
                   'avgnegvalue':avgvaluenegtrend,'numneg':numneg,'avgneg':avgneg, 'negpercent':avgnegpercent})
writer = pd.ExcelWriter('SpatialInfoTrendAveragess.xlsx', engine='xlsxwriter')  

df.to_excel(writer, sheet_name='Sheet1')
writer.save()  



#%% plot all maps

#want to use basemap here to also plot outline of states....
#states_file ='tl_2017_us_state/tl_2017_us_state'
#
#lat = dataset.variables['lat'][:]
#lon = dataset.variables['lon'][:]
#
##lon_alt = 360-lon[::-1]
#lon_alt = lon-360
#
#
#for a_ind,a in enumerate(allresults):
#    if a_ind<10:
#        
#        #fig = plt.figure(figsize=(8, 8))
#
#        m = Basemap(projection='merc',llcrnrlat=np.min(lat),urcrnrlat=np.max(lat),\
#            llcrnrlon=np.min(lon_alt),urcrnrlon=np.max(lon_alt))
#        
#        
#        
#        m.readshapefile(states_file,'tl_2017_us_state',drawbounds=True,color='0.3')
#        #x,y = np.meshgrid(lon_alt,lat)
#        
#        lons, lats = m.makegrid(300, 120) # get lat/lons of ny by nx evenly space grid.
#        x, y = m(lons, lats) 
#        
#        
#        cmap = plt.cm.gist_rainbow
#        #cmap.set_under ('1.0')
#        #cmap.set_bad('0.8')
#        
##        m.pcolormesh(x,y,a,cmap='RdBu_r')
##        plt.clim(-.1, .1)
#        
#        ext =[np.min(x),np.max(x),np.min(y),np.max(y)]
#        
#        a[a==0]=np.nan
#                
#        fname = filenames[a_ind]
#        
#        if 'corr' in fname:
#            m.imshow(a,extent = ext,cmap='bwr')
#            #plt.clim(-.1,.1)
#        else:            
#            m.imshow(a,extent = ext,cmap='Reds')
#            #plt.clim(0,.1)
#
#
#        # plot the states shape file
#
#        
##        
##        fig = plt.figure(a_ind,frameon=False)
##        ax = fig.add_axes([0, 0, 1, 1])
##        ax.axis('off')
##        a[a == 0.0] = np.nan
##        ax.imshow(a,origin='lower',vmin = -1, vmax=1, cmap = 'RdBu')
#        plt.title(filenames[a_ind])
##        plt.show()
        
    
    

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