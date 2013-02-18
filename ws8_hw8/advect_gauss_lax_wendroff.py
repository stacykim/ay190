import sys,math
from pylab import *


# THE LAX-WENDROFF ADVECTION EQUATION SOLVER -------------------------------------

def lax_wendroff(y):
    """
    Implements the Lax-Wendroff method to solve the advection equation.  The
    solution is not computed at the boundary points (call apply_bcs to do this).
    """
    return [y[i]-v*dt/(2*dx)*(y[i+1]-y[i-1])+0.5*(v*dt/dx)**2*(y[i-1]-2*y[i]+y[i+1])
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

        y=lax_wendroff(y)    # integrate advection equation
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
                imin=min(-20+v*t,len(x))
                imax=min(25+v*t,len(x))
                xlim([imin,imax])
                ylim([-0.005,0.01])
                #yscale('log')
                text(10,0.8*max(max(y),max(yana)),'t={0}'.format(t))

                draw()    
        except:
            print t
       """ 
        if t%150==0 and t>=150 and t<=600:
            nplot+=3
            subplot(4,3,nplot)
            plot(x,y,'b-')     # plot numerical result
            plot(x,yana,'r-')  # plot analytic results
            ha = 60 if (nplot <=3 or idt==2) else 10
            va2 = 0.4 if idt==2 else 0.6
            text(ha,0.8*max(max(y),max(yana)),'t={0}'.format(t))
            text(ha,va2*max(max(y),max(yana)),'cfl={0}'.format(cfl[idt]))
            draw()


#savefig('lax_wendroff_snapshots.pdf')
show()
ioff()
clf()

for idt in range(len(cfl)-1):
    dt=cfl[idt]*dx/abs(v)
    t=arange(0,dt*(int)(tf/dt+1),dt)
    plot(t,err[idt])

xlabel('t')
ylabel('error')
ylim([0,0.0002])
legend(['cfl={0}'.format(cfl[0]),'cfl={0}'.format(cfl[1]),
        'cfl={0}'.format(cfl[2])],loc='best')
savefig('lax_wendroff_err.pdf')
show()

err2=[err[0][i] for i in range(len(err[0])) if (i%2==0 and not i==1)]
dt=cfl[1]
t=arange(0,dt*(int)(tf/dt+1),dt)
print len(t), len(err[1]), len(err2)
plot(t,array(err[1])/array(err2))
savefig('lax_wendroff_convg.pdf')
show()
