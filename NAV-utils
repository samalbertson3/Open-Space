from math import *
import csv

def lng_rot(lng):
    return [[cos(lng),-sin(lng),0],[sin(lng),cos(lng),0],[0,0,1]]

def lat_rot(lat):
    return [[1,0,0],[0,cos(-lat),-sin(-lat)],[0,sin(-lat),cos(-lat)]] #Had to switch signs. Why?

def calculate(r0,lng0,lat,v0,theta0,phi0,t):
    #Extracts orbital elements from r,v vectors
    GM = 3.5316000*(10**12)
    r0 = r0 + 600000.
    theta = deg_to_rad(theta0)
    phi = deg_to_rad(phi0)
    #lng = lng0
    lng = adj_lng(lng0,t)
    r = [r0*sin(lng)*cos(lat),-r0*cos(lng)*cos(lat),r0*sin(lat)]
    v = [v0*sin(theta)*cos(phi),-v0*sin(phi),v0*cos(theta)*cos(phi)]
    v = mat_mult(lat_rot(lat),v)
    v = mat_mult(lng_rot(lng),v)
    h = cross(r,v)
    n = cross([0,0,1],h)
    #e = mult(sub(mult(r,squ(v)-GM/r0),mult(v,dot(r,v))),1/GM)
    e = sub(mult(cross(v,h),1/GM), mult(r,1/r0))
    p = squ(h)/GM
    a = p/(1-abs_val(e)**2)
    i = acos(h[2]/abs_val(h))
    meridian = [0,-1,0]
    #meridian = mat_mult(lng_rot(lng),meridian)
    if n[0] < 0:
        omega = 2*pi-acos(dot(meridian,n)/(abs_val(n)))
    else:
        omega = acos(dot(meridian,n)/(abs_val(n)))
    w = acos(dot(n,e)/(abs_val(n)*abs_val(e)))
    nu = acos(dot(e,r)/(abs_val(e)*abs_val(r)))
    return a,e,rad_to_deg(i),rad_to_deg(omega),rad_to_deg(w),rad_to_deg(nu)

def iterate(start,step,stop):
    #Creates a number line vector from [start] to [stop] in step size [step]
    out = [start]
    if step > 0:
        while out[-1] < stop:
            out.append(out[-1]+step)
    if step < 0:
        while out[-1] > stop:
            out.append(out[-1]-step)
    return out



def adj_lng(lng,time):
    #for a ship at a given lng, find lng of that point at T=0
    out = (pi + lng + time/P*2*pi)%(2*pi) - pi
    #out = lng - 2*pi*((time/P)%P)
    return out

def conv_time(Y,D,h,m,s):
    return s + 60*m + 3600*h + D*P + 9203545*Y

def state_vectors(r0,lng0,lat,v0,theta0,phi0,t):
    r0 = r0 + 600000.
    theta = deg_to_rad(theta0)
    phi = deg_to_rad(phi0)
    lng = lng0
    #lng = adj_lng(lng0,t)
    r = [r0*sin(lng)*cos(lat),-r0*cos(lng)*cos(lat),r0*sin(lat)]
    v = [v0*sin(theta)*cos(phi),-v0*sin(phi),v0*cos(theta)*cos(phi)]
    lat_rot = [[1,0,0],[0,cos(-lat),-sin(-lat)],[0,sin(-lat),cos(-lat)]] #Had to switch signs. Why?
    lng_rot = [[cos(lng),-sin(lng),0],[sin(lng),cos(lng),0],[0,0,1]]
    v = mat_mult(lat_rot,v)
    v = mat_mult(lng_rot,v)
    return r,v

##def heading(a,e,i,omega,w,nu):
##    e = abs_val(e)
##    i = deg_to_rad(i)
##    omega = deg_to_rad(omega)
##    w = deg_to_rad(w)
##    nu = deg_to_rad(nu)
##    GM = 3.5316000*(10**12)
##    r = a*(1-e**2)/(1+e*cos(nu)) - 600000.
##    r_prime = sqrt(GM*a)/(a*(1-e**2))*e*sin(nu)
##    rt_prime = sqrt(GM*a)*(1+e*cos(nu))/(a*(1-e**2))
##    v = sqrt((r_prime)**2+(rt_prime)**2)
##    phi = atan(r_prime/rt_prime)
##    theta = i*cos(nu+w)
##    lng = (omega + nu + w)
##    lat = (nu+w)*sin(i)
##    if lng > pi:
##        lng = -pi + (lng-pi)
##    if lat > pi:
##        lat = -pi + (lat-pi)
##    return [r,rad_to_deg(lng),rad_to_deg(lat),v,rad_to_deg(theta),rad_to_deg(phi)]

def heading(a,e,i,omega,w,nu,t):
    e = abs_val(e)
    i = deg_to_rad(i)
    omega = deg_to_rad(omega)
    w = deg_to_rad(w)
    nu = deg_to_rad(nu)

    #switched from Rx(i) to Ry(-i); why does this work?
    r_hat = mat_mult(Rz(omega),mat_mult(Ry(-i),mat_mult(Rz(nu+w),[0,-1,0])))
    r = a*(1-e**2)/(1+e*cos(nu))
    r = r - 600000.

    lat = acos(dot([r_hat[0],r_hat[1],0],r_hat)/(abs_val([r_hat[0],r_hat[1],0])*abs_val(r_hat)))
    if r_hat[2] < 0:
        lat = -lat
    lng = acos(dot([0,-1,0],[r_hat[0],r_hat[1],0])/(abs_val([r_hat[0],r_hat[1],0])))
    if r_hat[0] < 0:
        lng = -lng
    lng = adj_lng(lng,-t)
    

    theta_hat = mat_mult(Rz(omega),mat_mult(Ry(-i),mat_mult(Rz(nu+w),[1,0,0])))
    r_prime = sqrt(GM*a)/(a*(1-e**2))*e*sin(nu)
    rt_prime = sqrt(GM*a)*(1+e*cos(nu))/(a*(1-e**2))
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

##def state_vectors2(a,e,i,omega,w,nu):
##    i = deg_to_rad(i)
##    omega = deg_to_rad(omega)
##    w = deg_to_rad(w)
##    nu = deg_to_rad(nu)
##    e = abs_val(e)
##    r0 = [cos(i)*sin(nu+w),-cos(i)*cos(w+nu),sin(i)*sin(w+nu)]
##    #M = [[cos(omega),-sin(omega),0],[sin(omega),cos(omega),0],[0,0,1]]
##    M = [[1,0,0],[0,1,0],[0,0,1]]
##    r = mult(mat_mult(M,r0),a*(1-e**2)/(1+e*cos(nu)))
##    return r
    

ap_ref = 600000.+589431.
pe_ref = 600000.+94065.
a_ref = (ap_ref+pe_ref)/2.
e_ref = ap_ref/a_ref-1
p_ref = a_ref*(1-e_ref**2)

r_err = 1
lng_err = 0.01
lat_err = 0.01
v_err = 0.1
theta_err = 1
phi_err = 5

#Kerbin stats
P = 5*3600 + 59*60 + 9.4
GM = 3.5316000*(10**12)
#386339,-2.917,0.605,1846.9,113,15,1088727
#545765,-2.865,0.25 ,1553.9,129,7,1089156
#589432,-2.612,-0.051,1479.3,131,0,1089522
#462085,-2.07,-0.548,1703,118,-15,1090151
#155690,-0.598,-0.508,2365.7,60,-10,1090784
#94131,0.104,0.071,2534.9,49,0,1091058
#225702,1.266,0.702,2191.8,82,13,1091463
#381068,2.006,0.613,1857.3,112,15,1091765

