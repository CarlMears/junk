#########################################
# This is from lusitanica/gauss-kruger-python-pt on github
#
# modified by Carl Mears to use numpy
######################################################


#from math import *
import numpy as np

############################################################################
# Meridian Arc
# Based on https://en.wikipedia.org/wiki/Meridian_arc#Expansions_in_the_eccentricity_(e)
# Expansions in the third flattening
############################################################################
def arcmer(a,equad,lat1,lat2):
    b=a*np.sqrt(1-equad)                #scalar 
    n=(a-b)/(a+b)                       #scalar - third flattening
    equad_test = 4.0*n/((1.0+n)*(1.0+n))
    print(equad,equad_test)
    H0=1.+((n**2)/4.)+((n**4)/64.)      #scalar
    H2=-(3./2.)*n + (3.0/16.0)*(n**3)   #scalar
    H4=(15./16.)*(n**2) - (15.0/64.0)*(n**4)   #scalar
    H6=-(35./48.)*(n**3)                 #scalar
    H8=(315.0/512.0)*(n**4)
    
    s1=((a+b)/2.0)*(H0*lat1 + H2*np.sin(2.*lat1) + H4*np.sin(4.*lat1)+H6*np.sin(6.*lat1))# + H8*np.sin(8.*lat1))
    s2=((a+b)/2.0)*(H0*lat2 + H2*np.sin(2.*lat2) + H4*np.sin(4.*lat2)+H6*np.sin(6.*lat2))# + H8*np.sin(8.*lat2))
    return s2-s1
#############################################################################
# Gauss-Kruger Projection
#############################################################################
def geogauss(lat,lon,a,equad,lat0,lon0):

    lat0=np.radians(lat0)
    lon0=np.radians(lon0)

    lat=np.radians(lat)
    lon=np.radians(lon)

    lon=lon-lon0

    N=a/np.sqrt(1-equad*(np.sin(lat))**2)
    
    RO=a*(1-equad)/((1-equad*(np.sin(lat)**2))**(3./2.))

    k1=(N/RO)+(4.*(N**2)/(RO**2))-((np.tan(lat))**2)

    k2=(N/RO)-((np.tan(lat))**2)

    k3=N/RO*(14.-58.*((np.tan(lat)))**2)+40.*((np.tan(lat))**2)+((np.tan(lat))**4)-9.

    x=lon*N*np.cos(lat)+(lon**3)/6.*N*((np.cos(lat))**3)*k2+(lon**5)/120.*N*((np.cos(lat))**5)*k3

    y=arcmer(a,equad,lat0,lat)+(lon**2)/2.*N*np.sin(lat)*np.cos(lat)+((lon**4)/24.)*N*np.sin(lat)*((np.cos(lat))**3)*k1

    return x,y

#############################################################################
# Local Gauss Kruger
#############################################################################
def local_geogauss(lat0,lon0,lat,lon):
    
    # GRS-84
    a = 6378137.0  #semi major axis
    equad =0.00669437999  #first eccentricity squared 	
    coords = geogauss(lat,lon,a,equad,lat0,lon0)
    x,y = coords
    return x,y
#############################################################################
# PT-TM06/ETRS-89
#############################################################################
# def pttm(lat,lon):
#     # GRS-80
#     a = 6378137.  #semi major axix
#     equad =0.00669437999  #first eccentricity squared
#     # Natural Origin 
#     lat0=(39.+(40./60.)+(5.73/3600.))
#     lon0=(-(8.+(7./60.)+(59.19/3600.))) 	
#     coords = geogauss(lat,lon,a,equad,lat0,lon0)
#     x,y = coords
#     return x,y
#############################################################################
# Test values 
#############################################################################    

if __name__ == "__main__":

    a = 6378137.0  #semi major axis
    equad =0.00669437999  #first eccentricity squared 
    lat1 = np.radians(-0.5)
    lat2 = np.radians(0.5)
    dy_equator = arcmer(a,equad,lat1,lat2)
    print(f"1 deg lat at equator = {dy_equator}")

    lat0 = 0.0
    lon0 = 0.0

    dlat = 1.0
    dlon = 1.0

    lat = np.array([lat0,lat0,lat0+dlat],dtype=np.float64)
    lon = np.array([lon0,lon0+dlon,lon0],dtype=np.float64)

    x,y = local_geogauss(lat0,lon0,lat,lon)

    print(f"1 deg lon at equator = {(x[1]-x[0])/1000.0}")
    print(f"1 deg lon from eq length = {a*2.0*np.pi/360.0}")
    print(f"1 deg lat at equator = {(y[2]-y[0])/1000.0}")
    print(lon)
    print(lat)
    print(x)
    print(y)

