from math import *

def deg_to_rad(deg):
    rad = pi/180.*deg
    return rad
def rad_to_deg(rad):
    deg = 180./pi*rad
    return deg
def cross(r,v):
    out = [0,0,0]
    out[0] = r[1]*v[2]-r[2]*v[1]
    out[1] = r[2]*v[0]-r[0]*v[2]
    out[2] = r[0]*v[1]-r[1]*v[0]
    return out
def squ(x):
    tot = 0
    for i in x:
        tot += i**2
    return tot
def mult(x,y):
    z = range(0,len(x))
    for i in range(0,len(x)):
        z[i] = x[i]*y
    return z
def dot(x,y):
    z = range(0,len(x))
    for i in range(0,len(x)):
        z[i] = x[i]*y[i]
    return sum(z)
def sub(x,y):
    z = range(0,len(x))
    for i in range(0,len(x)):
        z[i] = x[i]-y[i]
    return z
def abs_val(x):
    tot = 0
    for i in x:
        tot += i**2
    return sqrt(tot)
def mat_mult(M,v):
    out = [0,0,0]
    out[0] = M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2]
    out[1] = M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2]
    out[2] = M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2]
    return out

def mean(mat):
    return sum(mat)/len(mat)

