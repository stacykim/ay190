# set10.py
# created 2/20/13 by stacy kim

from numpy import *
from nbodyBH import *
from matplotlib.pyplot import *
import sys

iter=0

def gravi_nbody(r,p1,Npm1):
    """
    Equations of motion for particle p1 under the gravitational forces from
    the other Npm1 particles.  Calculates only the velocity component.
    """
    #x=array([x,y,z])

    xsum,ysum,zsum=0,0,0
    for p2 in Npm1:
        dist=(r[0]-p2.r[0])**2+(r[1]-p2.r[1])**2+(r[2]-p2.r[2])**2+a**2
        xsum+=p2.m*(r[0]-p2.r[0])/dist**1.5
        ysum+=p2.m*(r[1]-p2.r[1])/dist**1.5
        zsum+=p2.m*(r[2]-p2.r[2])/dist**1.5
        
    return -G*array([xsum,ysum,zsum])


# MAIN ---------------------------------------------------------------------------

a=0         # force-softening constant
ntest=100   # number of particles to compare direct and BH methods

fn='HW5_Data.txt'
N=len(open(fn,'rU').read().split('\n')[1:])-2


# Initialize particles and oct-tree
particles=init_particles_from_file(fn,N)
root=create_root(particles)
root.subdivide()
print 'done subdividing'


# directly calculate force on particles
fg_direct=[]
start_time=time()
for i in range(ntest):
    p=particles[i]
    if i<5: print p.r
    fg_direct.append(gravi_nbody(p.r,p,particles[:i]+particles[i+1:]))
time_direct=time()-start_time
print 'took',time_direct,'s to compute forces directly'

np_fg_direct=array(fg_direct)
for i in range(5): print np_fg_direct[i],norm(np_fg_direct[i])


# calculate force via BH algorithm using various opening angles
opening_angles=[0.01,0.02,0.05,0.075,0.1,0.125,0.15,0.2,0.25,0.3]
noa=len(opening_angles)
times=zeros(noa)
err=zeros(noa)
for j in range(noa):
    opening_angle=opening_angles[j]
    start_time=time()
    fg_BH=[]

    for i in range(ntest):
        p=particles[i]
        fg_BH.append(gravi_nbodyBH(p.r,p,root,opening_angle,a))

    times[j]=time()-start_time
    np_fg_BH=array(fg_BH)
    for i in range(5): print np_fg_BH[i],norm(np_fg_BH[i])
    err[j]=sum([norm(np_fg_BH[i])/norm(np_fg_direct[i]) for i in range(ntest)])/ntest

    print 'oa={0}: {1} sec with average error of {2}'.format(opening_angle,times[j],err[j])


# algorithm performance comparisons
plot(opening_angles,times)
plot([0,0.3],[time_direct,time_direct])
xlabel('opening angle')
ylabel('runtime (s)')
savefig('runtimes.pdf')
show()

plot(opening_angles,err)
xlabel('opening angle')
ylabel('average accuracy')
savefig('accuracy.pdf')
show()

plot(opening_angles,[1-ee for ee in err])
xlabel('opening angle')
ylabel('average error')
savefig('errors.pdf')
show()

