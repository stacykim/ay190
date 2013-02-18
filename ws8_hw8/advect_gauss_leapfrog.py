import sys,math
from pylab import *


# THE LEAPFROG ADVECTION EQUATION SOLVER ----------------------------------------

def leapfrog(y,yold2):
    """
    Implements the leapfrog method to solve the advection equation.  The solution
    is not computed at the boundary points (call apply_bcs to do this).
    """
    return [yold2[i] - v*dt/dx*(y[i+1]-y[i-1]) for i in range(1,len(y)-1)]

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
    nplot=0

    for it in range(ntmax):
        t = it*dt
        # save previous and previous previous data
        yold2 = yold
        yold = y

        y=leapfrog(y,yold2)   # integrate advection equation
        y=apply_bcs(x,y)     # after update, apply boundary conditions
        yana = analytic(x,t) # get analytic result for time t
        
        # compute error estimate
        err[idt][it] = sqrt(sum([(yana[i]-y[i])**2 for i in range(len(x))]))/n

        """
        try:
            if t%5==0 and t>=5:
                clf()              # clear figure and ...
                plot(x,y,'b-')     # plot numerical result
                plot(x,yana,'r-')  # plot analytic results
                #imin=min(-20+v*t,len(x))
                #imax=min(25+v*t,len(x))
                #xlim([imin,imax])
                #ylim([-0.005,0.01])
                #yscale('log')
                text(10,0.8*max(max(y),max(yana)),'t={0}'.format(t))

                draw()    
        except:
            print t
        """
        if idt==0 and ((t%30==0 and t>=30 and t<=240) or \
                       ((t+10)%30==0 and t>=700 and t<=800)):
            nplot+=1
            subplot(3,4,nplot)
            plot(x,y,'b-')     # plot numerical result
            plot(x,yana,'r-')  # plot analytic results
            imin=min(-20+v*t,len(x))
            imax=min(25+v*t,len(x))
            xlim([imin,imax])
            ylim([-0.005,0.01])
            text(imin+5,0.006,'t={0}'.format(t))
            draw()


savefig('leapfrog_snapshots.pdf')
show()
ioff()
clf()
colors = ['b','g','r']
for idt in range(len(cfl)):
    dt=cfl[idt]*dx/abs(v)
    t=arange(0,dt*(int)(tf/dt+1),dt)
    plot(t,err[idt])

xlabel('t')
ylabel('error')
ylim([0,0.001])
legend(['cfl={0}'.format(cfl[0]),'cfl={0}'.format(cfl[1]),
        'cfl={0}'.format(cfl[2])],loc='best')
savefig('leapfrog_err.pdf')
show()
