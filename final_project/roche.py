# roche.py
# created by stacy kim

from numpy import *
from findroot import *
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D


G = 6.6733e-8     # gravitational constant
MSUN = 1.9891e33  # mass of the sun (g)
AU = 1.49598e13
DAYS = 3600*24   # seconds in a day

err=2e26

M1 = 40*MSUN # mass of primary
M2 = 20*MSUN # mass of secondary

sep=0.2*AU         # distance between primary and secondary
a = M2*sep/(M1+M2) # location of primary
b = M1*sep/(M1+M2) # location of secondary
omega = 2*pi / (5.5998 * DAYS) # angular velocity of orbiting bodies

V = lambda x2, y2: -(G*M1/sqrt((sqrt(x2)-a)**2 + y2)
                    + G*M2/sqrt((sqrt(x2)+b)**2 + y2)
                    + omega**2 * (x2 + y2)/2.)


# =============================================================================
# PLOT POTENTIAL
fig = figure()
ax = fig.add_subplot(111,projection='3d')
xlabel('AU')
ylabel('AU')

v=[]
mi=-0.75
ma=0.75
st=0.02
for x in arange(mi*AU,ma*AU,st*AU):
    for y in arange(mi*AU,ma*AU,st*AU):
        v+=[V(x*x,y*y)]
n=len(arange(mi,ma,st))
y=concatenate([arange(mi,ma,st)]*n)
x=concatenate([ones(n)*y[i] for i in range(n)])
ax.plot(x,y,v,'ko',ms=1)
savefig('potential.pdf')
show()


# =============================================================================
# PLOT ROCHE SURFACES

def f(x,y,v):
    r1 = sqrt( (x-a)**2 + y**2 )
    r2 = sqrt( (x+b)**2 + y**2 )
    x3 = omega**2 * r1 * r2 * (x**2 + y**2)/2.
    return r1*r2*v + G*(M1*r2 + M2*r1) + x3

mid=-0.0488
err = 2e26
surfaces=[]
for v in [-5e15+1e4,-5e15,-5.22e15,-7.5e15,-1e16]:
    print v
    st = 0.005 #if v != -5.22e15 else 0.0001
    if v <= -5*10**15: mi,ma = -0.25,0.25
    else:              mi,ma = -0.60,0.60

    # compute surface for negative x
    xsurface=[]
    ysurface=[]
    for x in arange(mi*AU,mid*AU,st*AU):
        ff = lambda y: f(x,y,v)
        try:
            yy,cnt=bisection(0,ma*AU,ff,err)
            xsurface+=[x]
            ysurface+=[yy]
        except ValueError:  pass

    surfaces+=[[array(xsurface)/AU,array(ysurface)/AU]]

    # calculate surface for positive x separately
    xsurface=[]
    ysurface=[]
    for x in arange(mid*AU,ma*AU,st*AU):
        ff = lambda y: f(x,y,v)
        try:
            yy,cnt=bisection(0,ma*AU,ff,err)
            xsurface+=[x]
            ysurface+=[yy]
        except ValueError:  pass

    surfaces+=[[array(xsurface)/AU,array(ysurface)/AU]]


# plot 2D surfaces
for i,surface in enumerate(surfaces):    
    plot(surface[0], surface[1],'k')
    plot(surface[0],-surface[1],'k')

    # connect edges of surfaces across y=0
    if not (i%2==1 and i<=3):
        plot( [surface[0][0]]*2, [surface[1][0],-surface[1][0]], 'k')
    if not (i%2==0 and i<=3):
        plot( [surface[0][-1]]*2, [surface[1][-1],-surface[1][-1]], 'k')

ylim([-0.4,0.4])
savefig('roche_surfaces2D.pdf')
show()


# plot 3D surfaces
fig=figure()
ax = fig.add_subplot(111,projection='3d')
n=len(arange(0,pi,pi/20))

yz = concatenate([surfaces[0][1],-surfaces[0][1]])
z = concatenate([yz*sin(th) for th in arange(0,pi,pi/20)])
y = concatenate([yz*cos(th) for th in arange(0,pi,pi/20)])
x  = concatenate([surfaces[0][0]]*n*2)
plot(x,y,z,'ko',ms=1)

yz = concatenate([surfaces[1][1],-surfaces[1][1]])
z = concatenate([yz*sin(th) for th in arange(0,pi,pi/20)])
y = concatenate([yz*cos(th) for th in arange(0,pi,pi/20)])
x  = concatenate([surfaces[1][0]]*n*2)
plot(x,y,z,'ko',ms=1)

ylim(xlim())
ax.set_zlim(xlim())

show()
