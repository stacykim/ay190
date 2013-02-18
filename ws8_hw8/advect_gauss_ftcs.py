import sys,math
from pylab import *


# THE FTCS ADVECTION EQUATION SOLVER --------------------------------------------

def ftcs(y):
    """
    Implements the FTCS method to solve the advection equation.  The solution is
    not computed at the boundary points (call apply_bcs to do this).
    """
    return [y[i]-v*dt/(2*dx)*(y[i+1]-y[i-1]) for i in range(1,len(y)-1)]

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
cfl = array([1.0])
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
    nplot=0

    for it in range(ntmax):
        t = it*dt
        # save previous and previous previous data
        yold2 = yold
        yold = y

        y=ftcs(y)        # integrate advection equation
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
        if cfl[idt]==1.0 and t%25==0 and t>=25 and t<=225:
            nplot+=1
            subplot(3,3,nplot)
            plot(x,y,'b-')     # plot numerical result
            plot(x,yana,'r-')  # plot analytic results
            pos= 'right' #if nplot <= 6 else 'left'
            text(60,0.9*max(y),'t={0}'.format(t))
            draw()

savefig('ftcs_snapshots.pdf')
ioff()
show()
clf()
idt=0
dt=cfl[idt]*dx/abs(v)
t=arange(0,dt*(int)(tf/dt+1),dt)
plot(t[:103],err[idt][:103])
xlabel('t')
ylabel('error')
savefig('ftcs_err.pdf')
show()

for i in range(200): print '{0} {1}'.format(t[i],err[idt][i])
