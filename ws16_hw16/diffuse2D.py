#!/usr/bin/env python
#from __future__ import division
from matplotlib.patches import Patch
from pylab import *


# ================================================================================
# SUPPORTING FUNCTIONS

def gauss(x,y,x0,y0,sig):
    """
    Returns a 2D gaussian centered at (x0, y0) and with width sig.
    Assumes that x and y are in the form returned by numpy.meshgrid().
    Used to set initial condition.
    """
    return exp( -((x-x0)**2 + (y-y0)**2)/(2*sig**2) ) / (sig*sqrt(2*pi))


def ring(x,y):
    """
    Returns a ring of width 0.05 stretching from 0.05 to 0.1 from the center at
    (0.5,0.5).  Assumes that x and y are in the form returned by numpy.meshgrid().
    Used to set intiial condition.
    """
    val=sqrt( (x-0.5)**2 + (y-0.5)**2 )
    z=array([array([1.0 if 0.05 <= r <= 0.1 else 0.0 for r in row]) for row in val])
    return z


def calc_ddxy(Z,nx,ny,dx2):
    """
    Returns the spatial second derivative of the diffusion equation, using the
    centered finite difference method.
    """
    ddy=array([(Z[i,2:]-2*Z[i,1:-1]+Z[i,:-2])/(dx*dx) for i in range(1,nx-1)])
    ddx=array([(Z[2:,i]-2*Z[1:-1,i]+Z[:-2,i])/(dx*dx) for i in range(1,ny-1)])

    return ddy+ddx.T


def set_bound(Z,nx,ny):
    # boundaries
    Z[0,0] = Z[1,1]
    Z[0,:] = Z[1,:]
    Z[:,0] = Z[:,1]
    Z[nx-1,ny-1] = Z[nx-2,nx-2]
    Z[nx-1,:] = Z[nx-2,:]
    Z[:,ny-1] = Z[:,ny-2]

    return Z


# ================================================================================
# SOLVE THE 2D DIFFUSION EQUATION

D = 1.0
nit = 1000
nx = 100
ny = 100
"""
x0 = 50
y0 = 50
sig = 10.0
x = linspace(0, 100, nx)
y = linspace(0, 100, ny)
X,Y = meshgrid(x, y)
Z = gauss(X,Y,x0,y0,sig)
"""
x = linspace(0,1,nx)
y = linspace(0,1,nx)
X,Y = meshgrid(x,y)
Z = ring(X,Y)

dx = x[1]-x[0]
dx2 = dx*dx
idx2 = 1.0/dx2
dt = dx2/(4*D) # max timestep before FTCS becomes unstable
print 'dx =',dx,'\ndt =',dt

RHS = zeros((nx,ny))

ion()
pcolor(X, Y, Z, vmin=0.0, vmax=0.5)
colorbar()
nplot=0
for it in range(nit):
    Z[1:-1,1:-1] += D*calc_ddxy(Z,nx,ny,dx2)*dt
    Z = set_bound(Z,nx,ny)

    """
    if it % 50 == 0:
        clf()
        pcolor(X, Y, Z, vmin=0.0, vmax=0.5)
        colorbar()
        text(10,80,'t={0}'.format(round(it*dt,3)))
        #fname = 'frame%04d.png'%it
        #print 'Saving frame', fname
        #savefig(fname)
        draw()
    """
    if it % 90 == 0 and nplot < 15:
        nplot+= 1 if (nplot+1) % 4 != 0 else 2
        subplot(4,4,nplot,aspect='equal')
        pcolor(X, Y, Z, vmin=0.0, vmax=0.5)
        axis('off')
        text(0.1,0.8,'t={0}'.format(round(it*dt,3)))
        draw()

ioff()
savefig('2Dring_evol.pdf')
show()
