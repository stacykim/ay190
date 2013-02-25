# particles.py
# created 2/27/12 by stacykim (cut class from nbody.py)
#
# Implements a particle with mass, position, and velocity.
#
#
# modified 3/4/12: cast x, v as Vectors upon initialization
# modified 2/10/13: 
#     switched support of vector arith. from homegrown Vector class to numpy


from numpy import *


class Particle(object):
    """A particle with mass m, position r and velocity v."""

    def __init__(self, m, r, v):
        #r=[x,y,z], v=[vx,vy,vz]
        self.m=m
        self.r=array(r)
        self.v=array(v)
