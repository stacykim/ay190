# interp.py
# created 1/14/13 by stacy kim
#
# Contains a suite of interpolation routines.
#
# NOTE: In all of the folowing routines, x is an array of values to interpolate the
#  discretized function given by arrays xdat, ydat.
#
# modified 3/7/13: pw_quadratic and pw_linear returned interpolated values in a 
#     2D array of single-element arrays; should now return 1D array

import sys, math
import numpy as np
from operator import mul


def lagrange(xdat, ydat, x):
    """Lagrange interpolation."""
    xdat=np.array(xdat,dtype=np.float)
    ydat=np.array(ydat,dtype=np.float)

    return [sum([ydat[j]*reduce(mul,[(xx-xk)/(xdat[j]-xk) for xk in np.concatenate((xdat[:j],xdat[j+1:]))])
                 for j in range(len(xdat))]) for xx in x]


def pw_linear(xdat, ydat, x):
    """Piecewise linear interpolation (for non-piecewise, supply 2 data pts)."""
    # sort data (x,y) pairs by xdat values
    a=np.array([xdat,ydat],dtype=np.float)
    dat=a.T[a.T[:,0].argsort()]
    xsrt,ysrt=dat[:,0],dat[:,1]

    # interpolate at given locations x
    y=[] # interpolated values
    for xx in x:
        for i in range(len(xsrt)):
            if xx < xsrt[i]: break
        i-=1
        y.append(lagrange(xsrt[i:i+2],ysrt[i:i+2],[xx])[0])
        #y.append(ysrt[i]+(ysrt[i+1]-ysrt[i])/(xsrt[i+1]-xsrt[i])*(xx-xsrt[i]))

    return y


def pw_quadratic(xdat, ydat, x):
    """Piecewise quadratic interpolation (for non-piecewise, supply 3 data pts)."""
    # sort data (x,y) pairs by xdat values
    a=np.array([xdat,ydat],dtype=np.float)
    dat=a.T[a.T[:,0].argsort()]
    xsrt,ysrt=dat[:,0],dat[:,1]

    # interpolate at given locations x
    y=[] # interpolated values
    for xx in x:
        for i in range(len(xsrt)-1):
            if xx< xsrt[i]: break
        i-= 1
        y.append(lagrange(xsrt[i:i+3],ysrt[i:i+3],[xx])[0])

    return y


def cubic_hermite(xdat, ydat, ypdat, x):
    """Piecewise cubic Hermite interpolation."""
    psi0 = lambda z: 2*z**3-3*z**2+1
    psi1 = lambda z: z**3-2*z**2+z

    # assumes xdat,ydat,ypdat are already sorted by xdat values
    xsrt,ysrt,ypsrt=xdat,ydat,ypdat

    # sort data (x,y,yp) triples by xdat values
    a=np.array([xdat,ydat,ypdat])
    dat=a.T[a.T[:,0].argsort()]
    xsrt,ysrt,ypsrt=dat[:,0],dat[:,1],dat[:,2]

    # interpolate at given locations x
    y=[] # interpolated values
    for xx in x:
        for i in range(len(xsrt)-1):
            if xx < xsrt[i]: break
        i-=1
        z=(xx-xsrt[i])/(xsrt[i+1]-xsrt[i])
        y.append(ysrt[i]*psi0(z)+ysrt[i+1]*psi0(1-z)
                 +ypsrt[i]*(xsrt[i+1]-xsrt[i])*psi1(z)
                 -ypsrt[i+1]*(xsrt[i+1]-xsrt[i])*psi1(1-z))

    return y
