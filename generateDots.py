# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 06:52:06 2017

@author: mattdzugan
"""
import random
import time
import json
import numpy as np
import matplotlib.pyplot as plt
import shapely.geometry as shp
import math
import cartopy.feature as cfeature
from cartopy.io import shapereader

thisRepo = '/home/mattdzugan/Documents/dev/dribbble debut/'
land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])

###############################################################################
###############################################################################
##                          Useful Functions                                 ##
###############################################################################
###############################################################################
def greatCircleDistance(lon1,lat1,lon2,lat2):
    deg2rad = math.pi/180.0
    R = 6378.1
    lon1r = lon1*deg2rad
    lat1r = lat1*deg2rad
    lon2r = lon2*deg2rad
    lat2r = lat2*deg2rad
    if (abs(lon1-lon2)+abs(lat1-lat2)>0.0000001):
        centralAngle = math.acos( math.sin(lat1r)*math.sin(lat2r) + math.cos(lat1r)*math.cos(lat2r)*math.cos(lon2r-lon1r)  )
    else:
        centralAngle = 0
    return centralAngle*R
    
def reckon(lon,lat,bearing,distance):
    R = 6378.1
    d = float(distance)
    b = math.radians(bearing)
    lat1 = math.radians(lat)
    lon1 = math.radians(lon)
    lat2 = math.asin( math.sin(lat1)*math.cos(d/R) +
                 math.cos(lat1)*math.sin(d/R)*math.cos(b))
    lon2 = lon1 + math.atan2(math.sin(b)*math.sin(d/R)*math.cos(lat1),
                 math.cos(d/R)-math.sin(lat1)*math.sin(lat2))
    lat2 = math.degrees(lat2)
    lon2 = math.degrees(lon2)    
    return [lon2,lat2]
    


def poissonDisc(points,R1,k):
    N = len(points)
    # determine how many points are active
    activePoints = [x for x in points if x['active'] == True]
    Na = len(activePoints)
    while Na>0:
        # there's still room to (try) adding more
        # pick a point
        nn = random.randint(0,Na-1)
        searchPoint = activePoints[nn]
        successYet = False
        # give it k attempts
        for ii in range(k):
            bearing = 360*random.random()
            distance = math.sqrt(3*R1*R1*random.random()+R1*R1)
            testlonlat = reckon(searchPoint['lon'],searchPoint['lat'],bearing,distance)
            tLon = testlonlat[0]
            tLat = testlonlat[1]
            tooClose = False
            for jj in range(N):
                if (abs(tLat-points[jj]['lat'])<(R1/100.0)): #only do rigorous check if they have a chance (really good speedup)
                    d = greatCircleDistance(tLon,tLat,points[jj]['lon'],points[jj]['lat'])
                    if (d<R1):
                        tooClose = True
                        break
            if not tooClose:
                successYet = True
                break
        # if no success mark the point inactive
        if successYet:
            newPoint = {'lat':tLat, 'lon':tLon, 'active':True}
            points.append(newPoint)
        else:
            activePoints[nn]['active'] = False
        # recalculate intro stuff
        N = len(points)
        # determine how many points are active
        activePoints = [x for x in points if x['active'] == True]
        Na = len(activePoints)
    return points


###############################################################################
###############################################################################
##                             Main Script                                   ##
###############################################################################
###############################################################################

## Definitions
pRad = 95.0 #radius of specs (in kilometers) (Use 95 for actual basketball scale)
nGon = 16    #the order polygon meant to approximate a circle

## Generate List of Points
points = []
firstPoint = {'lat':(-90+180.0*random.random()), 'lon':(-180+360.0*random.random()), 'active':True}
points.append(firstPoint)
start = time.time()
points = poissonDisc(points,pRad,30) 
end = time.time()
print(str(len(points))+" points in "+str(round(end-start))+" seconds")

## Clean up data
for point in points:
    point['land'] = False
    point['lon'] = point['lon']%360.0
    if point['lon']>180.0:
       point['lon'] -= 360.0 

## Determine if they're Land or Water
shpfilename = shapereader.natural_earth(resolution='110m',
                                        category='physical',
                                        name='land')
reader = shapereader.Reader(shpfilename)
countries = reader.records()
for country in countries:
    for ii in range(len(points)):
        if country.geometry.contains(shp.Point(points[ii]['lon'],points[ii]['lat'])):
            points[ii]['land'] = True

## For each point build a little circle
for ii in range(len(points)):
    pLat = points[ii]['lat']
    pLon = points[ii]['lon']
    vLats = []
    vLons = []
    #vP = []
    gJ = []
    for nn in range(nGon+1):
        brng = (360.0/nGon)*nn
        newLL = reckon(pLon,pLat,brng,pRad/2.0)
        vLat = newLL[1]
        vLon = newLL[0]
        vLon = vLon%360.0
        if vLon>180.0:
           vLon -= 360.0 
        vLon = round(vLon,2)
        vLat = round(vLat,2)
        vLons.append(vLon)
        vLats.append(vLat)
        #vP.append({'x':vLon, 'y':vLat})
        gJ.append([vLon, vLat])
    points[ii]['poly'] = [gJ]
    if points[ii]['land']:
        plt.plot(vLons, vLats,'g')
    else:
        plt.plot(vLons, vLats,'b')
plt.show()    


landList=np.array(list(o['land'] for o in points))
lonList=np.array(list(o['lon'] for o in points))
latList=np.array(list(o['lat'] for o in points))
plt.scatter(lonList, latList, s=100, c=landList)
plt.show()

## Build geojson
landGeoJson =  {"type":"Feature", "geometry":{"type":"MultiPolygon", "coordinates":[]}}
waterGeoJson = {"type":"Feature", "geometry":{"type":"MultiPolygon", "coordinates":[]}}
for point in points:
    if point['land']:
        landGeoJson['geometry']['coordinates'].append(point['poly'])
    else:
        waterGeoJson['geometry']['coordinates'].append(point['poly'])

featureCollection = {"type":"FeatureCollection", "features":[landGeoJson,waterGeoJson]}
## and export
with open(thisRepo+"data/dotPolys_"+str(int(pRad))+".geo.json", 'wb') as outfile:
        json.dump(featureCollection, outfile)