from math import *
import numpy

#Longitudinal rotation matrix
def lng_rot(lng):
    return numpy.matrix([[cos(lng),-sin(lng),0],[sin(lng),cos(lng),0],[0,0,1]])

#Latitudinal rotation matrix
def lat_rot(lat):
    return numpy.matrix([[1,0,0],[0,cos(-lat),-sin(-lat)],[0,sin(-lat),cos(-lat)]]) #Had to switch signs. Why?)

def calculate(r0,lng0,lat,v0,theta0,phi0,t):
    #Extracts Keplerian orbital elements from heading data
    #Assumes orbit around Kerbin
    #lng0, lat are input in rad
    #theta0, phi0 are input in deg
    #r0 is input in m
    #v0 is input in m/s
    GM = 3.5316000*(10**12)
    r0 = r0 + 600000.
    theta = deg_to_rad(theta0)
    phi = deg_to_rad(phi0)
    lng = adj_lng(lng0,t)
    r = numpy.matrix([r0*sin(lng)*cos(lat),-r0*cos(lng)*cos(lat),r0*sin(lat)])
    v = numpy.matrix([v0*sin(theta)*cos(phi),-v0*sin(phi),v0*cos(theta)*cos(phi)])
    v = lat_rot(lat)*numpy.transpose(v) #latitude, then longitude rotations
    v = lng_rot(lng)*v #order is important!
    h = numpy.cross(r,numpy.transpose(v))
    n = numpy.cross(numpy.matrix([0,00,1]),h)
    e = numpy.cross(numpy.transpose(v),h)/GM - r/numpy.linalg.norm(r)
    p = numpy.linalg.norm(h*numpy.transpose(h))/GM
    a = p/(1-numpy.linalg.norm(e)**2)
    i = acos(h[0,2]/numpy.linalg.norm(h))
    meridian = numpy.matrix([0,-1,0]) #Meridian points in direction of 0 degrees longitude
    if n[0,0] < 0:
        omega = 2*pi-acos(meridian*numpy.transpose(n)/(numpy.linalg.norm(n)))
    else:
        omega = acos(meridian*numpy.transpose(n)/(numpy.linalg.norm(n)))
    w = acos(n*numpy.transpose(e)/(numpy.linalg.norm(n)*numpy.linalg.norm(e)))
    nu = acos(e*numpy.transpose(r)/(numpy.linalg.norm(e)*numpy.linalg.norm(r)))
    return a,e,rad_to_deg(i),rad_to_deg(omega),rad_to_deg(w),rad_to_deg(nu)

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
    #Converts heading data to r,v state vectors in reference frame at epoch
    r0 = r0 + 600000.
    theta = deg_to_rad(theta0)
    phi = deg_to_rad(phi0)
    lng = adj_lng(lng0,t)
    r = numpy.matrix([r0*sin(lng)*cos(lat),-r0*cos(lng)*cos(lat),r0*sin(lat)])
    v = numpy.matrix([v0*sin(theta)*cos(phi),-v0*sin(phi),v0*cos(theta)*cos(phi)])
    v = lat_rot(lat)*numpy.transpose(v) #latitude, then longitude rotations
    v = lng_rot(lng)*v #order is important!
    return r,v

def heading(a,e,i,omega,w,nu,t):
    #Converts Keplerian orbital elements to heading
    i = deg_to_rad(i)
    omega = deg_to_rad(omega)
    w = deg_to_rad(w)
    nu = deg_to_rad(nu)
    m = numpy.transpose(numpy.matrix([0,-1,0]))
    r_hat = Rz(omega)*Ry(-i)*Rz(nu+w)*m
    r = a*(1-e**2)/(1+e*cos(nu))
    r = r - 600000.
    lat = acos(numpy.matrix([r_hat[0,0],r_hat[1,0],0])*r_hat/(numpy.linalg.norm([r_hat[0,0],r_hat[1,0],0])*numpy.linalg.norm(r_hat)))
    if r_hat[2,0] < 0:
        lat = -lat
    lng = acos(numpy.transpose(m)*numpy.transpose(numpy.matrix([r_hat[0,0],r_hat[1,0],0]))/(numpy.linalg.norm([r_hat[0,0],r_hat[1,0],0])))
    if r_hat[0,0] < 0:
        lng = -lng
    lng = adj_lng(lng,-t)
    

    theta_hat = Rz(omega)*Ry(-i)*Rz(nu+w)*numpy.transpose(numpy.matrix([1,0,0]))
    r_prime = sqrt(GM*a)/(a*(1-e**2))*e*sin(nu)
    rt_prime = sqrt(GM*a)*(1+e*cos(nu))/(a*(1-e**2))
    print 'r\' = ', r_prime
    print 'rv\' = ', rt_prime
    v = sqrt(r_prime**2 + rt_prime**2) ####

    theta = atan(theta_hat[0,0]/theta_hat[1,0]) #####
    phi = atan(e*sin(nu)/(1+e*cos(nu)))
    return r,lng,lat,v,rad_to_deg(theta),rad_to_deg(phi)

def Rx(t):
    #Rotate around x-axis
    return numpy.matrix([[1,0,0],[0,cos(t),-sin(t)],[0,sin(t),cos(t)]])

def Ry(t):
    #Rotate around y-axis
    return numpy.matrix([[cos(t),0,sin(t)],[0,1,0],[-sin(t),0,cos(t)]])

def Rz(t):
    #Rotate around z-axis
    return numpy.matrix([[cos(t),-sin(t),0],[sin(t),cos(t),0],[0,0,1]])

def deg_to_rad(deg):
    rad = pi/180.*deg
    return rad
def rad_to_deg(rad):
    deg = 180./pi*rad
    return deg

#Measurement errors in heading
r_err = 1 #m
lng_err = 0.01 #rad
lat_err = 0.01 #rad
v_err = 0.1 #m/s
theta_err = 1 #deg
phi_err = 5 #deg

#Kerbin stats
P = 5*3600 + 59*60 + 9.4
GM = 3.5316000*(10**12)

#Sample heading data
#386339,-2.917,0.605,1846.9,113,15,1088727
#545765,-2.865,0.25 ,1553.9,129,7,1089156
#589432,-2.612,-0.051,1479.3,131,0,1089522
#462085,-2.07,-0.548,1703,118,-15,1090151
#155690,-0.598,-0.508,2365.7,60,-10,1090784
#94131,0.104,0.071,2534.9,49,0,1091058
#225702,1.266,0.702,2191.8,82,13,1091463
#381068,2.006,0.613,1857.3,112,15,1091765
