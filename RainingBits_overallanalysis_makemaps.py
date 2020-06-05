#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Updated May 2020

This code analyzes results from overall (70 year) analysis from code
CPC_overall_ITanalysis.py

inputs:
results file spatial_results_overall.pickle from CPC_overall_ITanalysis.py code
file of 'overall_Hx_values.pickle' from CPC_HXonly.py code, only contains precipitation entropies
a single precipitation file, e.g. precip.V1.0.1948.nc from CPC dataset
threshold values in Spatial_thresholds_mm.pickle computed from preivous code

outputs: all results saved in two file types
results_directional_averages.pickle                  
SpatialInfo_Averages.nc4

@author: allisongoodwell
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pickle

#from mpl_toolkits.basemap import Basemap

precip_folder='C:\\Users\\goodwela\\Dropbox\\UCDenver\\research\\rainfall_research\\CPC_raingage_gridded_dailydata\\original_precip_data\\'
string = precip_folder + 'precip.V1.0.1948.nc'
dataset = Dataset(string)

np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

maxdim = 9 #9 lagged neighbors (including self) and 1 current target at center

lags = [1]
N = 2 #number of bins (2 for binary pdf)
years = range(1948,2018)
nyrs = len(years)
    
pptdata = np.asarray(dataset.variables['precip'][:])
lat = dataset.variables['lat'][:]
lon = dataset.variables['lon'][:]
t = dataset.variables['time'][:]
    
#vector resolution: find x and y components of directions
#angles = np.asfarray([0, 45, 90, 135, 180, 225, 270, 315])
#change to angle from north (clockwise)
angles = np.asfarray([90, 135, 180, 225, 270, 315, 0, 45])
angles = angles*np.pi/180

# open files:
result_file = 'spatial_results_overall.pickle'
f_myfile = open(result_file,'rb')
AllResults = pickle.load(f_myfile)
Long_Lists = pickle.load(f_myfile)
Lat = pickle.load(f_myfile)
f_myfile.close()

Hx_file = 'overall_Hx_values.pickle'
f_myfile = open(Hx_file,'rb')
Hxvals = pickle.load(f_myfile)
f_myfile.close()  

thresh_file='Spatial_thresholds_mm.pickle'                  
f_myfile = open(thresh_file,'rb')
threshvals = pickle.load(f_myfile)
f_myfile.close() 
   
                
#%% 
    
I_2dvals = np.empty((np.size(lat),np.size(lon),4,9))*np.nan
I_3dvals = np.empty((np.size(lat),np.size(lon),4,9))*np.nan
Eta_vals = np.empty((np.size(lat),np.size(lon),4,9))*np.nan
I_mainangle = np.empty((np.size(lat),np.size(lon),4))*np.nan
I_vectorstrength = np.empty((np.size(lat),np.size(lon),4))*np.nan
Eta_mainangle = np.empty((np.size(lat),np.size(lon),4))*np.nan
Eta_vectorstrength = np.empty((np.size(lat),np.size(lon),4))*np.nan

Eta_max_mainangle = np.empty((np.size(lat),np.size(lon),4))*np.nan
Eta_max_vectorstrength = np.empty((np.size(lat),np.size(lon),4))*np.nan

I_max_mainangle = np.empty((np.size(lat),np.size(lon),4))*np.nan
I_max_vectorstrength = np.empty((np.size(lat),np.size(lon),4))*np.nan

I_top4sources = np.empty((np.size(lat),np.size(lon),4))*np.nan

U_EW= np.empty((np.size(lat),np.size(lon),4))*np.nan
U_NS= np.empty((np.size(lat),np.size(lon),4))*np.nan
U_NESW= np.empty((np.size(lat),np.size(lon),4))*np.nan
U_NWSE= np.empty((np.size(lat),np.size(lon),4))*np.nan

S_EW= np.empty((np.size(lat),np.size(lon),4))*np.nan
S_NS= np.empty((np.size(lat),np.size(lon),4))*np.nan
S_NESW= np.empty((np.size(lat),np.size(lon),4))*np.nan
S_NWSE= np.empty((np.size(lat),np.size(lon),4))*np.nan

R_EW= np.empty((np.size(lat),np.size(lon),4))*np.nan
R_NS= np.empty((np.size(lat),np.size(lon),4))*np.nan
R_NESW= np.empty((np.size(lat),np.size(lon),4))*np.nan
R_NWSE= np.empty((np.size(lat),np.size(lon),4))*np.nan

I_fraction= np.empty((np.size(lat),np.size(lon),4))*np.nan
Eta_fraction= np.empty((np.size(lat),np.size(lon),4))*np.nan
I4_fraction= np.empty((np.size(lat),np.size(lon),4))*np.nan
I2_NS_fraction= np.empty((np.size(lat),np.size(lon),4))*np.nan
I2_EW_fraction= np.empty((np.size(lat),np.size(lon),4))*np.nan
        
#now need to do some filling to obtain maps from lists of dictionaries with results
#each item in AllResults is associated with a latitude from Lat
#then each item in AllResults[0] is associated with a longitude from  Long_lists[0]
#each item in AllResults[0][0] is a season 
    #AllResults[0][0][0]['I_2d'] should access a specific measure
lats = np.asfarray(lat)
lons = np.asfarray(lon)

I_test = np.empty((np.size(lat),np.size(lon)))*np.nan
    
for a_ind,A in enumerate(AllResults): #A is list of longitudes
    
    lat_value = Lat[a_ind]
    la_ind = np.where(lats==lat_value)
       
    for b_ind,B in enumerate(A): #B is single longitude
#                
        lon_value = Long_Lists[a_ind][b_ind]        
        lo_ind = np.where(lons==lon_value)
 

        for s, info in enumerate(B): #info is results structure for single season s
            
             
            #print('lat = '+ str(lat_value)+ ' lon = '+ str(lon_value)+ ' s=', str(s))
            
            I_2dvals[la_ind,lo_ind,s,:]=info['I_2d']
            I_3dvals[la_ind,lo_ind,s,:]=info['I_3d']
            Eta_vals[la_ind,lo_ind,s,:]=info['Eta']
            
            #here, want to omit points without fewer than 6 neighbors
            if np.sum(I_2dvals[la_ind,lo_ind,s,:]==0)>2:
                #print('skipping values'+str(la_ind)+' '+str(lo_ind))
                continue
            
            
            #compute primary direction (angle) and strength of info to center based on 2d metrics
            Ivals = info['I_2d']
            Ivals = Ivals[0:8] #index 8 is lagged middle cell itself, 9 is current value of middle cell
            
            
            I_yval = np.sum([np.cos(angles[i])*Ivals[i] for i in range(0,8)])
            I_xval = np.sum([np.sin(angles[i])*Ivals[i] for i in range(0,8)])
            
            I_vectorstrength[la_ind,lo_ind,s]=np.sqrt(I_xval**2 + I_yval**2)
            if I_xval>0:
                I_mainangle[la_ind,lo_ind,s] = np.arctan2(I_xval,I_yval)
            else:
                I_mainangle[la_ind,lo_ind,s] = 2*np.pi + np.arctan2(I_xval,I_yval)
            
            I_mainangle[la_ind,lo_ind,s]=I_mainangle[la_ind,lo_ind,s]*180/np.pi
            
            #find max I and asscoiated angle (without vector resolution)
            I_maxval = np.max(Ivals)
            I_maxindex = np.argmax(Ivals)
            I_angle = angles[I_maxindex]*180/np.pi
            I_max_mainangle[la_ind,lo_ind,s]=I_angle
            I_max_vectorstrength[la_ind,lo_ind,s]=I_maxval
            
            #print('Ivect strength = '+str(I_vectorstrength[la_ind,lo_ind,s]))
            #compute primary direction (angle) and strength of info to center based on TE metrics
            Etavals = info['Eta']
            Etavals = Etavals[0:8] #index 8 is lagged middle cell itself, 9 is current value of middle cell
            
            I_xval = np.sum([np.cos(angles[i])*Etavals[i] for i in range(0,8)])
            I_yval = np.sum([np.sin(angles[i])*Etavals[i] for i in range(0,8)])
            
            Eta_vectorstrength[la_ind,lo_ind,s]=np.sqrt(I_xval**2 + I_yval**2)
            if I_xval>0:
                Eta_mainangle[la_ind,lo_ind,s] = np.arctan2(I_xval,I_yval)
            else:
                Eta_mainangle[la_ind,lo_ind,s] = 2*np.pi + np.arctan2(I_xval,I_yval)
            
            Eta_mainangle[la_ind,lo_ind,s]=Eta_mainangle[la_ind,lo_ind,s]*180/np.pi
            
            
            
            #find max Eta and asscoiated angle (without vector resolution)
            Eta_maxval = np.max(Etavals)
            Eta_maxindex = np.argmax(Etavals)
            Eta_angle = angles[Eta_maxindex]*180/np.pi
            Eta_max_mainangle[la_ind,lo_ind,s]=Eta_angle
            Eta_max_vectorstrength[la_ind,lo_ind,s]=Eta_maxval
            
                       
            #next retain the totol 4-d info from top 4 U+S sources
            I_top4sources[la_ind,lo_ind,s]=info['I_top4']
            
            #next retain U1,U2,R,S, and Itotal from opposite direction sources
            
            #each has 4 different pieces...E-W, NE-SW, N-S, NW-SE
            I_2directions = np.asarray(info['I_directions'])
            U_sources = (np.asarray(info['UN'])+np.asarray(info['US']))/I_2directions
            R_sources = np.asarray(info['R'])/I_2directions
            S_sources = np.asarray(info['S'])/I_2directions
                
            U_EW[la_ind,lo_ind,s]=U_sources[0]
            U_NESW[la_ind,lo_ind,s]=U_sources[1]
            U_NS[la_ind,lo_ind,s]=U_sources[2]
            U_NWSE[la_ind,lo_ind,s]=U_sources[3]
            
            S_EW[la_ind,lo_ind,s]=S_sources[0]
            S_NESW[la_ind,lo_ind,s]=S_sources[1]
            S_NS[la_ind,lo_ind,s]=S_sources[2]
            S_NWSE[la_ind,lo_ind,s]=S_sources[3]
            
            R_EW[la_ind,lo_ind,s]=R_sources[0]
            R_NESW[la_ind,lo_ind,s]=R_sources[1]
            R_NS[la_ind,lo_ind,s]=R_sources[2]
            R_NWSE[la_ind,lo_ind,s]=R_sources[3]
                               
            I_fraction[la_ind,lo_ind,s]= I_vectorstrength[la_ind,lo_ind,s]/Hxvals[la_ind,lo_ind,s]
            Eta_fraction[la_ind,lo_ind,s]= Eta_vectorstrength[la_ind,lo_ind,s]/Hxvals[la_ind,lo_ind,s]
            I4_fraction[la_ind,lo_ind,s]= I_top4sources[la_ind,lo_ind,s]/Hxvals[la_ind,lo_ind,s]
            I2_NS_fraction[la_ind,lo_ind,s] = I_2directions[2]/Hxvals[la_ind,lo_ind,s]
            I2_EW_fraction[la_ind,lo_ind,s] = I_2directions[0]/Hxvals[la_ind,lo_ind,s]
               
print('done with information theory analysis....')
            
#%% loop over seasons, then info measure, then lat/lon to save maps
                                        
measures = [I_vectorstrength, I_mainangle, Eta_vectorstrength, Eta_mainangle, Eta_max_vectorstrength, Eta_max_mainangle,
            I_max_vectorstrength, I_max_mainangle,
            I_top4sources, U_EW, U_NESW, U_NS, U_NWSE,S_EW,S_NESW,S_NS, S_NWSE,
            R_EW, R_NESW, R_NS, R_NWSE, I_fraction, Eta_fraction, I4_fraction,I2_EW_fraction,I2_NS_fraction] 
measure_string =  ['Ivectorstrength', 'Imainangle', 'Etavectorstrength', 'Etamainangle',  'Eta_max_vectorstrength', 'Eta_max_mainangle',
                   'I_max_vectorstrength','I_max_mainangle',
            'Itop4', 'U_EW', 'U_NESW', 'U_NS', 'U_NWSE','S_EW','S_NESW','S_NS', 'S_NWSE',
            'R_EW', 'R_NESW', 'R_NS', 'R_NWSE','I_fraction','Eta_fraction','I4_fraction','I2_EW_fraction','I2_NS_fraction']

season_string = ['_DJF','_MAM','_JJA','_SON'] 


filenames=[]
allresults=[]

for m_ind,m in enumerate(measures):
    
    print('now doing analysis for: '+ measure_string[m_ind])
    
    for s in range(0,4):
        
        data = m[:,:,s]
        
        filenames.append(measure_string[m_ind]+season_string[s]+'_avg')
        allresults.append(data)



#%%
#plot histograms of fractions (I_fraction, Eta_fraction, I4_fraction, I2_EW_fraction, I2_NS_fraction)
        #season by season, first have to remove nan values
        
season_string =['winter','spring','summer','fall']
for s in range(0,4):
    

    
    oldMI = I_fraction[:,:,s]
    nanvals = np.isnan(oldMI)
    notnanvals = ~nanvals
    newMI = oldMI[notnanvals]
    

    
    oldI2EW = I2_EW_fraction[:,:,s]
    nanvals = np.isnan(oldI2EW)
    notnanvals = ~nanvals
    newI2EW = oldI2EW[notnanvals]
    

    
    oldI2NS = I2_NS_fraction[:,:,s]
    nanvals = np.isnan(oldI2NS)
    notnanvals = ~nanvals
    newI2NS = oldI2NS[notnanvals]
    

    
    oldTE = Eta_fraction[:,:,s]
    nanvals = np.isnan(oldTE)
    notnanvals = ~nanvals
    newTE = oldTE[notnanvals]
    

    
    oldI4 = I4_fraction[:,:,s]
    nanvals = np.isnan(oldI4)
    notnanvals = ~nanvals
    newI4 = oldI4[notnanvals]

    
    plt.figure(1,figsize=(6,6))
    
    plt.subplot(2,2,s+1)
    plt.hist(newMI,bins = np.asarray(range(0,35,1))/100, alpha=0.5,density=True)
    plt.hist(newTE,bins = np.asarray(range(0,35,1))/100, alpha=0.5,density=True)
    plt.hist(newI2EW,bins = np.asarray(range(0,35,1))/100, alpha=0.5,density=True)
    plt.hist(newI4, bins = np.asarray(range(0,35,1))/100, alpha=0.5,density=True)
    
    plt.legend(['MI','TE','I2','I4'])
    plt.title(season_string[s])

       
oldMI = I_fraction
nanvals = np.isnan(oldMI)
notnanvals = ~nanvals
newMI = oldMI[notnanvals]


oldI2EW = I2_EW_fraction
nanvals = np.isnan(oldI2EW)
notnanvals = ~nanvals
newI2EW = oldI2EW[notnanvals]


oldI2NS = I2_NS_fraction
nanvals = np.isnan(oldI2NS)
notnanvals = ~nanvals
newI2NS = oldI2NS[notnanvals]


oldTE = Eta_fraction
nanvals = np.isnan(oldTE)
notnanvals = ~nanvals
newTE = oldTE[notnanvals]


oldI4 = I4_fraction
nanvals = np.isnan(oldI4)
notnanvals = ~nanvals
newI4 = oldI4[notnanvals]


plt.figure(2,figsize=(3.25,2))
plt.hist(newMI,bins = np.asarray(range(0,35,1))/100, alpha=0.5,density=True)
plt.hist(newTE,bins = np.asarray(range(0,35,1))/100, alpha=0.5,density=True)
plt.hist(newI2EW,bins = np.asarray(range(0,35,1))/100, alpha=0.5,density=True)
plt.hist(newI4, bins = np.asarray(range(0,35,1))/100, alpha=0.5,density=True)

plt.legend(['$MI$','$TE$','$I_{EW}$','$I_4$'],fontsize=12)
plt.xlabel('fraction of uncertainty reduced')
plt.ylabel('% of grid cells')
plt.savefig('Fig_measurePDFS.pdf',bbox_inches='tight')

#%%
plt.figure(3)
plt.imshow(I_fraction[:,:,1]-Eta_fraction[:,:,1])
plt.colorbar()
                                
#%% now have a big list of filenames and allresults (list of lat-lon matrices)

#save in a pickle file

site_filenamesave='results_directional_averages.pickle'                  
f_myfile = open(site_filenamesave,'wb')
pickle.dump(allresults, f_myfile)
pickle.dump(filenames, f_myfile)
pickle.dump(threshvals, f_myfile)
f_myfile.close()   


    
#%% save output also as netcdf file...
dsout = Dataset("SpatialInfo_Averages.nc4", "w", format="NETCDF4_CLASSIC")

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
for a_ind,a in enumerate(allresults):    
    X = dsout.createVariable(filenames[a_ind], np.float64, ('lat','lon',))  
    X[:,:] = a        

dsout.close()
     

#%% now want to look at "wind rose style plot" of I and Eta strength and direction

#import windrose
#from matplotlib import cm
#
#seasonname=['DJF','MAM','JJA','SON']
#
#fig=plt.figure()
#for i in range(1,5):
#    speeds = np.reshape(I_vectorstrength[:,:,i-1],(np.size(I_vectorstrength[:,:,i-1]),))
#    directions = np.reshape(I_mainangle[:,:,i-1],(np.size(I_mainangle[:,:,i-1]),))
#    
#    directions =[d for d_ind,d in enumerate(directions) if speeds[d_ind]>0 ]
#    #now in degrees from east
#    
#    speeds = [s for s in speeds if s>0]
#        
##    rect=[0.5+(i-1)*0.5,0.5,0.4,0.4] 
##    wa=windrose.WindroseAxes(fig, rect)
##    fig.add_axes(wa)
##    wa.bar(directions, speeds, normed=True, opening=0.8, edgecolor='white',bins=[.02, .04, .06, .08, .1, 1])
##    
##    #ax = fig.add_subplot(2, 2, i)
#    ax = windrose.WindroseAxes.from_ax()
#    #ax.bar(directions, speeds, normed=True, opening=0.8, edgecolor='white',bins=[.02, .04, .06, .08, .1, 1])
#    ax.contourf(directions, speeds, bins=np.arange(0,.15,0.01), cmap=cm.viridis, lw=3,nsector=32)
#    #ax.set_legend()
#    ax.grid(False)
#    plt.savefig('figs/I_directions_'+seasonname[i-1]+'.pdf')
#    
#fig=plt.figure()
#for i in range(1,5):
#    speeds = np.reshape(Eta_vectorstrength[:,:,i-1],(np.size(Eta_vectorstrength[:,:,i-1]),))
#    directions = np.reshape(Eta_mainangle[:,:,i-1],(np.size(Eta_mainangle[:,:,i-1]),))
#    
#    directions =[d for d_ind,d in enumerate(directions) if speeds[d_ind]>0 ]
#    #now in degrees from east
#    
#    speeds = [s for s in speeds if s>0]
#    ax = windrose.WindroseAxes.from_ax()
#    #ax.bar(directions, speeds, normed=True, opening=0.8, edgecolor='white',bins=[.01, .02, .03, .04, 1])
#    ax.contourf(directions, speeds, bins=np.arange(0,.1,.005), cmap=cm.viridis, lw=3,nsector=32)
#    #ax.set_legend()
#    ax.grid(False)
#    plt.savefig('figs/Eta_directions_'+seasonname[i-1]+'.pdf')  
#    fig=plt.figure()
#    
#for i in range(1,5):
#    speeds = np.reshape(Eta_max_vectorstrength[:,:,i-1],(np.size(Eta_max_vectorstrength[:,:,i-1]),))
#    directions = np.reshape(Eta_max_mainangle[:,:,i-1],(np.size(Eta_max_mainangle[:,:,i-1]),))
#    
#    directions =[d for d_ind,d in enumerate(directions) if speeds[d_ind]>0 ]
#    #now in degrees from east
#    
#    speeds = [s for s in speeds if s>0]
#    ax = windrose.WindroseAxes.from_ax()
#    #ax.bar(directions, speeds, normed=True, opening=0.8, edgecolor='white',bins=[.01, .02, .03, .04, 1])
#    ax.contourf(directions, speeds, bins=np.arange(0,.2,.1), cmap=cm.hot, lw=3,nsector=32)
#    #ax.set_legend()
#    ax.grid(False)
#    plt.savefig('figs/EtaMax_directions_'+seasonname[i-1]+'.pdf')  
#    
#for i in range(1,5):
#    speeds = np.reshape(I_max_vectorstrength[:,:,i-1],(np.size(I_max_vectorstrength[:,:,i-1]),))
#    directions = np.reshape(I_max_mainangle[:,:,i-1],(np.size(I_max_mainangle[:,:,i-1]),))
#    
#    directions =[d for d_ind,d in enumerate(directions) if speeds[d_ind]>0 ]
#    #now in degrees from east
#    
#    speeds = [s for s in speeds if s>0]
#    ax = windrose.WindroseAxes.from_ax()
#    #ax.bar(directions, speeds, normed=True, opening=0.8, edgecolor='white',bins=[.01, .02, .03, .04, 1])
#    ax.contourf(directions, speeds, bins=np.arange(0,.2,.1), cmap=cm.hot, lw=3,nsector=32)
#    #ax.set_legend()
#    ax.grid(False)
#    plt.savefig('figs/IMax_directions_'+seasonname[i-1]+'.pdf') 
#    
#    
##%% I also want spatial wind roses....for each HUC2 basin?
#    
##load HUC2 shapefile - take all data from each HUC2 and make windrose plot, save figures
#
#import shapefile
#from shapely.geometry import Polygon
#from shapely.geometry import Point
#
#altlon = 360-lon
#
#sf_HUC = shapefile.Reader('C:\\Users\\goodwela\\Dropbox\\UCDenver\\rainfall_research\\CPC_raingage_gridded_dailydata\\shapefiles\\All_HUC2_shapefile\\All_HUC2_shapefile')
#nshapes_H = len(sf_HUC.shapes())
#records_H = sf_HUC.records()
#
#HUCpolys = []
#HUCnames=[]
#HUCnums=[]
#site_HUCnum = np.zeros((len(lat),len(lon)))
#site_HUCname = np.empty((len(lat),len(lon)), dtype=object)
#
#HUC_I_mainangle_watersheds=[]
#HUC_I_vectorstrength_watersheds=[]
#HUC_Eta_mainangle_watersheds=[]
#HUC_Eta_vectorstrength_watersheds=[] 
#    
#for s in range(0,nshapes_H):
#    
#    HUC_I_mainangle=[]
#    HUC_I_vectorstrength=[]
#    HUC_Eta_mainangle=[]
#    HUC_Eta_vectorstrength=[]
#    
#    HUC_Longitude =[]
#    
#    shape_ex = sf_HUC.shape(s)
#    place_string = records_H[s][11]
#    place_num=records_H[s][10]
#    
#    x_lon = np.zeros((len(shape_ex.points),))
#    y_lat = np.zeros((len(shape_ex.points),))
#    coords = np.zeros((len(shape_ex.points),2))
#    for ip in range(len(shape_ex.points)):
#        x_lon[ip] = shape_ex.points[ip][0]
#        y_lat[ip] = shape_ex.points[ip][1]
#        coords[ip,0]=shape_ex.points[ip][0]
#        coords[ip,1]=shape_ex.points[ip][1]
#
#    #create polygon of that watershed shape
#    poly = Polygon((coords))
#    #poly = MultiPoint(coords).convex_hull
#    HUCpolys.append(poly)
#    HUCnames.append(place_string)
#    HUCnums.append(place_num)
#    
#    #only search in box near x_lon and y_lat points
#    lat_short = lat[lat>=np.min(y_lat)]
#    lat_short = lat_short[lat_short<=np.max(y_lat)]
#    
#    lon_short = altlon[altlon>=np.min(-x_lon)]
#    lon_short = lon_short[lon_short<=np.max(-x_lon)]
#    
#    for la_ind,la in enumerate(lat):
#        print('lat= '+str(la))
#        
#        
#        #long_list = Long_Lists[la_ind]
#        
#        for lo_ind,lo in enumerate(altlon):
#     
#            inside = HUCpolys[s].contains(Point((-lo,la)))
#            #print(inside, HUCnames[s])
#            
#            if inside == True: #and len(long_list)>0:
#                
#                
#                site_HUCnum[la_ind,lo_ind]=int(HUCnums[s])
#                site_HUCname[la_ind,lo_ind]=HUCnames[s]
#                
#                HUC_I_mainangle.append(I_mainangle[la_ind,lo_ind,2])
#                HUC_I_vectorstrength.append(I_vectorstrength[la_ind,lo_ind,2])
#                HUC_Eta_mainangle.append(Eta_mainangle[la_ind,lo_ind,2])
#                HUC_Eta_vectorstrength.append(Eta_vectorstrength[la_ind,lo_ind,2])
#
#
#                print('lat= '+str(la)+' lon= ',str(lo)+ ' HUC location: ', HUCnames[s])
#                continue
#            
#    fig=plt.figure()
#    speeds = HUC_I_vectorstrength
#    directions = HUC_I_mainangle
#    
#    directions =[d for d_ind,d in enumerate(directions) if speeds[d_ind]>0 ]
#    #now in degrees from east    
#    speeds = [s for s in speeds if s>0]
#        
#    ax = windrose.WindroseAxes.from_ax()
#    #ax.bar(directions, speeds, normed=True, opening=0.8, edgecolor='white',bins=[.02, .04, .06, .08, .1, 1])
#    ax.contourf(directions, speeds, bins=np.arange(0,.15,0.01), cmap=cm.viridis, lw=3,nsector=16)
#    #ax.set_legend()
#    ax.grid(False)
#    ax.set_yticks([])
#    ax.set_xticks([])
#    plt.savefig('figs/I_JJA_'+HUCnames[s]+'.pdf')
#    
#    fig=plt.figure()
#    speeds = HUC_Eta_vectorstrength
#    directions = HUC_Eta_mainangle
#    
#    directions =[d for d_ind,d in enumerate(directions) if speeds[d_ind]>0 ]
#    #now in degrees from east    
#    speeds = [s for s in speeds if s>0]
#        
#    ax = windrose.WindroseAxes.from_ax()
#    #ax.bar(directions, speeds, normed=True, opening=0.8, edgecolor='white',bins=[.02, .04, .06, .08, .1, 1])
#    ax.contourf(directions, speeds, bins=np.arange(0,.15,0.01), cmap=cm.viridis, lw=3,nsector=16)
#    #ax.set_legend()
#    ax.grid(False)
#    ax.set_yticks([])
#    ax.set_xticks([])
#    plt.savefig('figs/Eta_JJA_'+HUCnames[s]+'.pdf') 
#     
#    HUC_I_mainangle_watersheds.append(HUC_I_mainangle) 
#    HUC_I_vectorstrength_watersheds.append(HUC_I_vectorstrength) 
#    HUC_Eta_mainangle_watersheds.append(HUC_Eta_mainangle) 
#    HUC_Eta_vectorstrength_watersheds.append(HUC_Eta_vectorstrength) 
    


#%%
   




     