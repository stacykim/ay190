import sys,math
from pylab import *

import advect_gauss_upwind as upwind


# THE LAX-FRIEDRICH ADVECTION EQUATION SOLVER -----------------------------------

def lax_friedrich(y):
    """
    Implements the Lax-Friedrich method to solve the advection equation. The
    solution is not computed at the boundary points (call apply_bcs to do this).
    """
    return [0.5*(y[i+1]+y[i-1]) - v*dt/(2*dx)*(y[i+1]-y[i-1]) 
            for i in range(1,len(y)-1)]

def apply_bcs(x,y):
    """
    Applies 'outflow' boundary conditions. Assumes boundary points not
    included in y and adds them to y.
    """
    return [y[0]] + y + [y[-1]]

def analytic(x,t):
    """Analytic solution to advection equation for a Gaussian."""
    return exp(-(x-v*t-x0)**2/(2*sigma**2))


# MAIN ---------------------------------------------------------------------------

# parameters
dx = 0.1
v = 0.1
# set up grid
x = arange(0,100,dx)
n = len(x)
cfl = array([0.5,1.0,2.0])
t = 0.0
tf=800.0

# for initial data
sigma = sqrt(15.0)
x0 = 30.0

#set up initial conditions
y = analytic(x,0)
yold2 = y
yold = y

# evolve (and show evolution)
ion()  # "interaction on" - required for animation via draw()
figure()
ylim([0,1])
plot(x,y,'x-') # numerical data
plot(x,analytic(x,t),'r-') # analytic data
show()

err=[]
for idt in range(len(cfl)):
    #set up initial conditions
    y = analytic(x,0)
    yold2 = y
    yold = y

    dt = cfl[idt]*dx/abs(v)
    ntmax = (int)(tf/dt+1)
    err.append(zeros(ntmax))
    nplot=idt-2

    for it in range(ntmax):
        t = it*dt
        # save previous and previous previous data
        yold2 = yold
        yold = y

        y=lax_friedrich(y)   # integrate advection equation
        y=apply_bcs(x,y)     # after update, apply boundary conditions
        yana = analytic(x,t) # get analytic result for time t
        
        # compute error estimate
        err[idt][it] = sqrt(sum([(yana[i]-y[i])**2 for i in range(len(x))]))/n

        """
        try:
            if t<600 and t%10==0 and t>=10:
                clf()              # clear figure and ...
                plot(x,y,'b-')     # plot numerical result
                plot(x,yana,'r-')  # plot analytic results
                #yscale('log')
                draw()    
        except:
            print t
        """
        if t%150==0 and t>=150 and t<=600:
            nplot+=3
            subplot(4,3,nplot)
            plot(x,y,'b-')     # plot numerical result
            plot(x,yana,'r-')  # plot analytic results
            ha = 60 if nplot <=3 else 10
            va2 = 0.5 if idt==2 else 0.6
            text(ha,0.8*max(max(y),max(yana)),'t={0}'.format(t))
            text(ha,va2*max(max(y),max(yana)),'cfl={0}'.format(cfl[idt]))
            draw()


savefig('lax_friedrich_snapshots.pdf')
show()
ioff()
clf()
colors=['b','g','r']

for idt in range(len(cfl)):
    dt=cfl[idt]*dx/abs(v)
    t=arange(0,dt*(int)(tf/dt+1),dt)
    plot(t,err[idt],colors[idt])
ylim([0,0.0021])

for idt in range(len(cfl)):
    dt=cfl[idt]*dx/abs(v)
    t=arange(0,dt*(int)(tf/dt+1),dt)
    plot(t,upwind.err[idt],'{0}--'.format(colors[idt]))

xlabel('t')
ylabel('error')
ylim([0,0.0021])
legend(['cfl={0}'.format(cfl[0]),'cfl={0}'.format(cfl[1]),
        'cfl={0}'.format(cfl[2])],loc='best')
savefig('lax_friedrich_err.pdf')
show()
