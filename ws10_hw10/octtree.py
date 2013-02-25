# octtree.py
# created 2/20/13 by stacykim (modified 2D quadtree.py)
#
# Implements a 3D oct-tree containing particles.


import sys
import math
from numpy import *
from particles import *


class OctTree(object):
    """
    A rectangular oct-tree subdivided into cells each with 0 or 1 particles.
    To set up a subdivided oct-tree, call
        root.create_root(particles)
        root.subdivide()
    where particles is an array of Particles. The OctTree and its particles can be
    printed (to separate files) by calling
        root.print_tree().
    """

    def __init__(self, xrange, yrange, zrange, particles):
        """
        Initialize a cubical cell that covers the given range and contains the
        given particles.
        """
        self.xrange=array(xrange)
        self.yrange=array(yrange)
        self.zrange=array(zrange)
        self.size=xrange[1]-xrange[0]        # length of each side of cell
        self.m=sum([p.m for p in particles]) # total mass contained in the cell
        # calculate center of mass
        self.com=reduce(lambda com, p: com+p.m*p.r, particles,array([0,0,0]))/self.m if self.m!=0 else array([0,0,0])
        self.children=particles  # either a Particle or an OctTree


    def subdivide(self):
        """
        Divides the given oct-tree into 8 cells.  Assumes that oct.children 
        is an array of Particles, which is reassigned to an array of OctTrees.
        """
        if len(self.children)<2: return

        # Calculate ranges for each cell
        x0,x2=self.xrange
        x1=(x2+x0)/2.

        y0,y2=self.yrange
        y1=(y2+y0)/2.

        z0,z2=self.zrange
        z1=(z2+z0)/2.

        # Assign particles to appropriate cell
        particles_sub=[[],[],[],[],[],[],[],[]]

        for p in self.children:
            # cells are numbered in the same order as in the 2D Cartesian system,
            # first for [z1,z2] then [z0,z1]
            if   (x1 <= p.r[0] <= x2) and (y1 <= p.r[1] <= y2) and (z1 <= p.r[2] <= z2): particles_sub[0]+=[p]
            elif (x0 <= p.r[0] <= x1) and (y1 <= p.r[1] <= y2) and (z1 <= p.r[2] <= z2): particles_sub[1]+=[p]
            elif (x0 <= p.r[0] <= x1) and (y0 <= p.r[1] <= y1) and (z1 <= p.r[2] <= z2): particles_sub[2]+=[p]
            elif (x1 <= p.r[0] <= x2) and (y0 <= p.r[1] <= y1) and (z1 <= p.r[2] <= z2): particles_sub[3]+=[p]
            elif (x1 <= p.r[0] <= x2) and (y1 <= p.r[1] <= y2) and (z0 <= p.r[2] <= z1): particles_sub[4]+=[p]
            elif (x0 <= p.r[0] <= x1) and (y1 <= p.r[1] <= y2) and (z0 <= p.r[2] <= z1): particles_sub[5]+=[p]
            elif (x0 <= p.r[0] <= x1) and (y0 <= p.r[1] <= y1) and (z0 <= p.r[2] <= z1): particles_sub[6]+=[p]
            elif (x1 <= p.r[0] <= x2) and (y0 <= p.r[1] <= y1) and (z0 <= p.r[2] <= z1): particles_sub[7]+=[p]
            else: 
                print 'Found particle outside cell:\n p.r={0}\n xrange={1}\n yrange={2}\n zrange={3}\n' \
                      ''.format(p.r,[x0,x1,x2],[y0,y1,y2],[z0,z1,z2])
                print len(self.children)
                for i in range(20): print self.children[i].r
                sys.exit()

        # Create 8 sub-OctTrees
        self.children=[OctTree([x1,x2],[y1,y2],[z1,z2],particles_sub[0]),
                       OctTree([x0,x1],[y1,y2],[z1,z2],particles_sub[1]),
                       OctTree([x0,x1],[y0,y1],[z1,z2],particles_sub[2]),
                       OctTree([x1,x2],[y0,y1],[z1,z2],particles_sub[3]),
                       OctTree([x1,x2],[y1,y2],[z0,z1],particles_sub[4]),
                       OctTree([x0,x1],[y1,y2],[z0,z1],particles_sub[5]),
                       OctTree([x0,x1],[y0,y1],[z0,z1],particles_sub[6]),
                       OctTree([x1,x2],[y0,y1],[z0,z1],particles_sub[7])]

        # Recursively divide sub-OctTrees
        for child in self.children: child.subdivide()


    def print_tree(self, fqn, fpn):
        # Open oct-tree output file
        fq=open(fqn, 'w')
        fq.write(fqn+'\n')

        # Open particles output file
        fp=open(fpn, 'w')
        fp.write(fpn+'\n')

        self.print_oct(fq, fp)

        fq.close()
        fp.close()


    def print_oct(self, fq, fp):
        """
        Prints the edges of the given OctTree and its sub-OctTrees to the given
        output file and the particles in the OctTrees to the given particle output 
        file.
        """
        s=''
        for x in self.xrange:
            for y in self.yrange:
                for z in self.yrange:
                    s+='{0} {1} {2} '.format(x,y,z)
        fq.write(s)

        for oct in self.children:
            if type(oct)==OctTree:   oct.print_oct(fq,fp)
            elif type(oct)==Particle: fp.write('{0} {1} {2}\n'.format(oct.r[0],oct.r[1],oct.r[2]))


def create_root(particles):
    """
    Creates a root OctTree with spatial ranges large enough to contain the
    furthest particle.
    """
    # Find spatial range of particles
    range=0
    for p in particles:
        m=max(abs(p.r))
        if m > range: range=m

    range += 1e-4  # so no particles exactly on border

    return OctTree([-range, range], [-range, range], [-range,range], particles)
