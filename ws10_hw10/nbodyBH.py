# nbodyBH.py
# created 2/27/12 by stacykim (imported body of nobody.py)
#
# Computes a gravitational N-body simulation using the Barnes-Hut algorithm
# with force softening.  The position of each particle is calculated using 
# the sympletic Euler method.
#
#
# modified 2/27/12: exported Particle class to particles.py
# modified 3/4/12:  implemented Barnes-Hut algorithm
# modified 3/6/12
#    added output of total energy at every timestep
#    added function that calculates norm of a given vector, norm()
# modified 2/20/13: modified to only support BH algorithm (i.e. removed script)


import sys
import math
from numpy import *
from random import random
from time import time
from particles import *
from octtree import *

G=1               # gravitational constant

# BARNES-HUT ALGORITHM ----------------------------------------------------------

depth=0
def gravi_nbodyBH(r, p, oct, opening_angle,a):
    """
    Equations of motion for the Particle p at updated position r under 
    the gravitational forces from the particles in the OctTree oct.  
    Recursively calculates only the velocity component.
    """
    #r=array([x,y,z])
    global depth

    depth+=1
    # Check if particle inside given cell
    inside=((oct.xrange[0] <= p.r[0] <= oct.xrange[1]) and
            (oct.yrange[0] <= p.r[1] <= oct.yrange[1]) and
            (oct.zrange[0] <= p.r[2] <= oct.zrange[1]))

    # If no particle inside, then no force from this cell
    if len(oct.children) == 0:
        return array([0,0,0])

    # If one particle, then calculate exact 2-particle force
    elif len(oct.children) == 1:
        if not inside: return two_body(r,oct.children[0],a)
        else:          return array([0,0,0])

    # If more than one particle, det if approx using COM or not
    elif len(oct.children) > 1:
        # Calculate proper com (i.e. if particle in cell, recalculate com without it)
        if not inside: com=oct.com
        else:          com=oct.m*(oct.com-p.m*p.r)/(oct.m-p.m)#oct.com-p.m*p.r/oct.m

        d=math.pow((r[0]-com[0])**2+(r[1]-com[1])**2+(r[2]-com[2])**2,0.5)
        if oct.size/d < opening_angle:
            return two_body(r,Particle(oct.m,com,[0,0,0]),a)
        else:
            accel = array([0,0,0])
            for cell in oct.children:
                accel = accel + gravi_nbodyBH(r,p,cell,opening_angle,a)
            return accel


def two_body(r,p,a):
    """
    Returns the force-softened gravitational acceleration imparted by
    Particle p on another particle at position x.
    """
    dist=(r[0]-p.r[0])**2+(r[1]-p.r[1])**2+(r[2]-p.r[2])**2+a**2
    ax=-G*p.m*(r[0]-p.r[0])/dist**1.5
    ay=-G*p.m*(r[1]-p.r[1])/dist**1.5
    az=-G*p.m*(r[2]-p.r[2])/dist**1.5

    return array([ax,ay,az])


def norm(x):
    """Returns the norm of the given 3-vector."""
    return math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])


def init_particles_from_file(fn,N):
    """
    Initializes and returns an array of the first N particles from the file with 
    name fn.  The file can either specify mass, position, and velocity of each 
    particle (one per line) or just the positions (in which case m=1 and v=[0,0,0]).
    """
    particles=[[]]*N
    # Read inputs from input file if given
    f=open(fn,'rU')
    f.readline()  # Skip first 3 lines of input file
    f.readline()

    # Read in initial conditions for first N particles given in file
    for i in range(N):
        line=f.readline()[:-1].rsplit(' ')
        for j in range(line.count('')): line.remove('')

        if len(line)==7:
            params = array([float(line[j]) for j in range(7)])
        elif len(line)==3: # assume only positions, set m=1,v=[0,0,0]
            params = array([1.]+[float(line[j]) for j in range(3)]+[0.,0.,0.])
        else:
            print '{0}: Unrecognized line in input file:{1}'.format(i,line)
            sys.exit()

        particles[i]=Particle(params[0],params[1:4],params[4:])

    f.close()
    return particles

    
def init_particles_random(r_range,v_range,dim,N):
    """
    Initializes and returns an array of N particles each with a random position in 
    a circle (dim==2) or sphere (dim==3) of radius r_range centered on the origin 
    and with random velocities in the range [-v_range, v_range].  The array of 
    particles is outputted in the file partiles<N>.dat.
    """
    fn='particles{0}.dat'.format(N)
    f=open(fn,'w')
    f.write('{1}\n{0:>19}'.format('m',fn))
    for i in range(3): f.write('r{0}'.format(i).rjust(19))
    for i in range(3): f.write('v{0}'.format(i).rjust(19))
    f.write('\n')

    particles=[[]]*N
    for i in range(N):
        mass = m if m != 0 else random()

        r=rrange*random()
        if dim==2:
            th=2*math.pi*random()
            x=[r*math.cos(th),r*math.sin(th),0]
        if dim==3:
            th=math.pi*random()
            ph=2*math.pi*random()
            x=[r*math.sin(th)*math.cos(ph),r*math.sin(th)*math.sin(ph),r*math.cos(th)]

        v=[2*vrange*random()-vrange for j in range(3)]
        if dim==2: v[2]=0

        particles[i]=Particle(mass,array(x),array(v))

        f.write('{0!s:>19}'.format(m))
        for j in range(3): f.write('{0!s:>19}'.format(x[j]))
        for j in range(3): f.write('{0!s:>19}'.format(v[j]))
        f.write('\n')

    f.close()
    return particles


def update_particle(p,root,h):
    """
    Calculate gravitational force on particle p from the other particles in the 
    OctTree root using the symplectic Euler method and a time step of h. Returns
    the updated position and velocity of the particle.
    """
    r=p.r+h*p.v
    v=p.v+h*gravi_nbodyBH(r,p,root)
    return r,v


def update(N,particles,h):
    """
    Evolves the N particles forward by a time step of h using the symplectic Euler
    method to integrate forward in time and the Barnes-Hut algorithm.
    """
    particles0=particles[:] # make a copy of the unevolved system
    root=create_root(particles0)
    root.subdivide()

    for i in range(N):
        p=particles0[i]
        particles[i].r,particles[i].v=update_particle(p,root,h)
