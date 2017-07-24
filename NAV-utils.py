from math import *
import numpy as np
import csv

#Longitudinal rotation matrix
def lng_rot(lng):
    return [[cos(lng),-sin(lng),0],[sin(lng),cos(lng),0],[0,0,1]]

#Latitudinal rotation matrix
def lat_rot(lat):
    return [[1,0,0],[0,cos(-lat),-sin(-lat)],[0,sin(-lat),cos(-lat)]] #Had to switch signs. Why? 

def navball_orbitalElementTransform(r0,lng0,lat,v0,theta0,phi0,t):
    #Extracts Keplerian orbital elements from navball data
    #Assumes orbit around Kerbin
    #lng0, lat are input in rad
    #theta0, phi0 are input in deg
    #r0 is input in m
    #v0 is input in m/s
    r0 = r0 + radius_Kerbin
    theta = deg_to_rad(theta0)
    phi = deg_to_rad(phi0)
    lng = adj_lng(lng0,t)
    r = [r0*sin(lng)*cos(lat),-r0*cos(lng)*cos(lat),r0*sin(lat)]
    v = [v0*sin(theta)*cos(phi),-v0*sin(phi),v0*cos(theta)*cos(phi)]
    v = mat_mult(lat_rot(lat),v) #latitude, then longitude rotations
    v = mat_mult(lng_rot(lng),v) #order is important!
    h = cross(r,v)
    n = cross([0,0,1],h)
    e = sub(mult(cross(v,h),1/GM_kerbin), mult(r,1/r0))
    p = squ(h)/GM_kerbin
    a = p/(1-abs_val(e)**2)
    i = acos(h[2]/abs_val(h))
    meridian = [0,-1,0] #Meridian points in direction of 0 degrees longitude
    if n[0] < 0:
        omega = 2*pi-acos(dot(meridian,n)/(abs_val(n)))
    else:
        omega = acos(dot(meridian,n)/(abs_val(n)))
    w = acos(dot(n,e)/(abs_val(n)*abs_val(e)))
    nu = acos(dot(e,r)/(abs_val(e)*abs_val(r)))
    return a,e,rad_to_deg(i),rad_to_deg(omega),rad_to_deg(w),rad_to_deg(nu)
    
def orbitalElement_navballTransform(a,e,i,omega,w,nu):
    #should find r and v vectors based on data produced from navball_orbitalElement transform
    omega=np.radians(omega)
    w=np.radians(w)
    nu=np.radians(nu)       #assuming that nu would be entered in degrees as that is what its returned in from the inverse function of this one
    p=a*(1-np.abs(e)**2)
    P=np.array[p_i(omega,w,i),p_j(omega,w,i),p_k(w,i)]         #collect variables
    Q=np.array[q_i(omega,w,i),q_j(omega,w,i),q_k(w,i)]
    radius = p/(1+e*np.cos(nu))
    radius_vector = radius*np.cos(nu*P)+radius*np.sin(radius*Q)      #unsure if this is taking the correct trig functions of these equations
    velocity_vector = np.sqrt(GM_kerbin/p)*(-np.sin(nu*P)+(e+np.cos(nu))*Q)  

    
    return radius_vector,velocity_vector
    
    #equations for finding components of unit vectors in perifocal coordinate system
    #source: https://en.wikipedia.org/wiki/Perifocal_coordinate_system
def p_i(omega,w,i):
    #assuming rad
    return np.cos(omega)*np.cos(w)-np.sin(omega)*np.cos(i)*np.sin(w)

def p_j(omega,w,i):
    #assuming rad
    return np.sin(omega)*np.cos(w)+np.cos(omega)*np.cos(i)*np.sin(w)

def p_k(w,i):
    #assuming rad
    return np.sin(i)*np.sin(w)
    
def q_i(omega,w,i):
    #assuming rad
    return -np.cos(omega)*np.sin(w)-np.sin(omega)*np.cos(i)*np.cos(w)
    
def q_j(omega,w,i):
    #assuming rad
    return -np.sin(omega)*np.sin(w)+np.cos(omega)*np.cos(i)*np.cos(w)
    
def q_k(w,i):
    #assuming rad
    return np.sin(i)*np.cos(w)
    
def adj_lng(lng,time):
    #Establishes a Kerbocentric reference frame at epoch UT=0
    #for a ship at a given lng, find lng of that point at epoch
    out = (pi + lng + time/P*2*pi)%(2*pi) - pi
    #out = lng - 2*pi*((time/P)%P)
    return out

def conv_time(Y,D,h,m,s):
    #Changes Y/D/h/m/s time to seconds
    return s + 60*m + 3600*h + D*P + 9203545*Y

def state_vectors(r0,lng0,lat,v0,theta0,phi0,t):
    #Converts navball data to r,v state vectors in reference frame at epoch
    r0 = r0 + radius_Kerbin
    theta = deg_to_rad(theta0)
    phi = deg_to_rad(phi0)
    lng = adj_lng(lng0,t)
    r = [r0*sin(lng)*cos(lat),-r0*cos(lng)*cos(lat),r0*sin(lat)]
    v = [v0*sin(theta)*cos(phi),-v0*sin(phi),v0*cos(theta)*cos(phi)]
    lat_rot = [[1,0,0],[0,cos(-lat),-sin(-lat)],[0,sin(-lat),cos(-lat)]] #Had to switch signs. Why?
    lng_rot = [[cos(lng),-sin(lng),0],[sin(lng),cos(lng),0],[0,0,1]]
    v = mat_mult(lat_rot,v)
    v = mat_mult(lng_rot,v)
    return r,v

def heading(a,e,i,omega,w,nu,t):
    #Converts Keplerian orbital elements to navball data
    e = abs_val(e)
    i = deg_to_rad(i)
    omega = deg_to_rad(omega)
    w = deg_to_rad(w)
    nu = deg_to_rad(nu)

    #switched from Rx(i) to Ry(-i); why does this work?
    r_hat = mat_mult(Rz(omega),mat_mult(Ry(-i),mat_mult(Rz(nu+w),[0,-1,0])))
    r = a*(1-e**2)/(1+e*cos(nu))
    r = r - radius_Kerbin

    lat = acos(dot([r_hat[0],r_hat[1],0],r_hat)/(abs_val([r_hat[0],r_hat[1],0])*abs_val(r_hat)))
    if r_hat[2] < 0:
        lat = -lat
    lng = acos(dot([0,-1,0],[r_hat[0],r_hat[1],0])/(abs_val([r_hat[0],r_hat[1],0])))
    if r_hat[0] < 0:
        lng = -lng
    lng = adj_lng(lng,-t)
    

    theta_hat = mat_mult(Rz(omega),mat_mult(Ry(-i),mat_mult(Rz(nu+w),[1,0,0])))
    r_prime = sqrt(GM_kerbin*a)/(a*(1-e**2))*e*sin(nu)
    rt_prime = sqrt(GM_kerbin*a)*(1+e*cos(nu))/(a*(1-e**2))
    print 'r\' = ', r_prime
    print 'rv\' = ', rt_prime
    v = sqrt(r_prime**2 + rt_prime**2) ####

    theta = atan(theta_hat[0]/theta_hat[1]) #####
    phi = atan(e*sin(nu)/(1+e*cos(nu)))
    return r,lng,lat,v,rad_to_deg(theta),rad_to_deg(phi)

def Rx(t):
    return [[1,0,0],[0,cos(t),-sin(t)],[0,sin(t),cos(t)]]

def Ry(t):
    return [[cos(t),0,sin(t)],[0,1,0],[-sin(t),0,cos(t)]]

def Rz(t):
    return [[cos(t),-sin(t),0],[sin(t),cos(t),0],[0,0,1]]

#Measurement errors in navball readings
r_err = 1 #m
lng_err = 0.01 #rad
lat_err = 0.01 #rad
v_err = 0.1 #m/s
theta_err = 1 #deg
phi_err = 5 #deg

#Kerbin stats					#need to add Mun stats
P = 5*3600 + 59*60 + 9.4		#sidereal rotation period in seconds
GM_kerbin = 3.5316000*(10**12) 	#m^3/s^2, gravitational parameter of kerbin
radius_Kerbin = 600000.			#meters, radius of Kerbin

#Sample navball data
#386339,-2.917,0.605,1846.9,113,15,1088727
#545765,-2.865,0.25 ,1553.9,129,7,1089156
#589432,-2.612,-0.051,1479.3,131,0,1089522
#462085,-2.07,-0.548,1703,118,-15,1090151
#155690,-0.598,-0.508,2365.7,60,-10,1090784
#94131,0.104,0.071,2534.9,49,0,1091058
#225702,1.266,0.702,2191.8,82,13,1091463
#381068,2.006,0.613,1857.3,112,15,1091765

