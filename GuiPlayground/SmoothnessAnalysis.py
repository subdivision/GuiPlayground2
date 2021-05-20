import math as m
import numpy as np

#-----------------------------------------------------------------------------
def compute_slopes(pts, b_open, subd_iter):
    N = len(pts)
    NN = N-1 if b_open else N
    slopes = []
    xs = [4.*float(i)/(2.**subd_iter) for i in range(NN)]
    for i in range(NN):
        curr = pts[i]
        next = pts[(i+1)%N]
        der = next - curr
        der /= np.linalg.norm(der)
        v = der[1]/der[0]
        if np.abs(v) > 10.:
            v = 0.
        slopes.append(v)
    return xs, slopes

#-----------------------------------------------------------------------------
def compute_curvature(pts, b_open, subd_iter):
    N = len(pts)
    NN = N-1 if b_open else N
    crvs = []
    xs = [4.*float(i)/(2.**subd_iter) for i in range(NN)]
    for i in range(NN):
        prev = pts[(i-1+N)%N]
        curr = pts[i]
        next = pts[(i+1)%N]
        p,t, q,u, s,z = prev[0], prev[1], curr[0], curr[1], next[0], next[1]
        v = ((p*u-q*t)*z**2+(-p*u**2+q*t**2-p*q**2+p**2*q)*z\
            +s*t*u**2+(-s*t**2+p*s**2-p**2*s)*u+(q**2*s-q*s**2)*t)\
            /((q-p)*z+(p-s)*u+(s-q)*t)
        crvs.append(v)
    return xs, crvs

#============================ END OF FILE ====================================
