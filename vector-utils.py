from math import *

#Useful vector functions to use when your computer refuses to install numpy

def deg_to_rad(deg):
    #Converts degrees to radians
    rad = pi/180.*deg
    return rad

def rad_to_deg(rad):
    #Converts radians to degrees
    deg = 180./pi*rad
    return deg

def cross(r,v):
    #Cross product of a 3x3 matrix
    out = [0,0,0]
    out[0] = r[1]*v[2]-r[2]*v[1]
    out[1] = r[2]*v[0]-r[0]*v[2]
    out[2] = r[0]*v[1]-r[1]*v[0]
    return out

def squ(x):
    #Sum-squared of vector elements
    tot = 0
    for i in x:
        tot += i**2
    return tot

def mult(x,y):
    #Multiplies vector X by scalar Y
    z = range(0,len(x))
    for i in range(0,len(x)):
        z[i] = x[i]*y
    return z

def dot(x,y):
    #Dot product of vectors X and Y
    z = range(0,len(x))
    for i in range(0,len(x)):
        z[i] = x[i]*y[i]
    return sum(z)

def sub(x,y):
    #Subtracts vector Y from vector X
    z = range(0,len(x))
    for i in range(0,len(x)):
        z[i] = x[i]-y[i]
    return z

def abs_val(x):
    #Vector magnitude
    tot = 0
    for i in x:
        tot += i**2
    return sqrt(tot)

def mat_mult(M,v):
    #Left-handed multiplication of vector V by matrix M
    out = [0,0,0]
    out[0] = M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2]
    out[1] = M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2]
    out[2] = M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2]
    return out

def mean(vec):
    #Mean of vector elements
    return sum(vec)/len(vec)

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
