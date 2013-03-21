from pylab import *
import sys,math

# ================================================================================
# SUPPORTING FUNCTIONS

def analytic(t,t0,x,x0,D):
    return (t0/(t0+t))**0.5 * exp( -(x-x0)**2 / (4*D*(t0+t)) )


def calc_rhs(D,y,dx):
    """
    Returns the spatial portion of the diffusion equation, calculating the second
    derivative of y using the centered finite difference method.
    """
    ddy=(y[2:]-2*y[1:-1]+y[:-2])/(dx*dx)
    return D*ddy*dt


def set_bound(y):
    y[0] = y[1]
    y[n-1] = y[n-2]
    return y


# ================================================================================
# SOLVE THE 1D DIFFUSION EQUATION

n = 200
x = linspace(0,100,n)
dx = x[1]-x[0]
dx2 = dx*dx

x0 = 50.0
D = 1.0
t0 = 1.0
t = 0.0
nt = 3000

# initial conditions
y = exp( -(x-x0)**2 / (4.0*D*t0) )
dt = dx2/(2*D) # max timestep before FTCS becomes unstable
print 'dt =',dt,'dx =',dx

nplot=0
ion()
plot(x,y,"r-")
plot(x,analytic(t,t0,x,x0,D),"bx-")
show()

ie=len(x)/2 # peak of gaussian where error analysis done
err_peak=zeros(nt)
err_side=zeros(nt)
for it in range(nt):
    y[1:-1] += calc_rhs(D,y,dx)
    y = set_bound(y)
    t += dt

    yana=analytic(t,t0,x,x0,D)

    #if it % 10 == 0:
    if it == int(10**(nplot/2.4)) and nplot <= 8:
        nplot+=1
        subplot(3,3,nplot)
        clf()
        plot(x,y,'r')
        plot(x,yana,'b')
        ylim([0,1])
        text(55,0.8,'t={0}'.format(round(t,2)))
        draw()

    err_peak[it]=abs(yana[ie]-y[ie])
    err_side[it]=abs(yana[0]-y[0])

ioff()
plot(x,y,"r")
savefig('1D_evol.pdf')
show()

plot(arange(0,nt*dt,dt),err_peak)
plot(arange(0,nt*dt,dt),err_side)
legend(['peak','boundary (x=0)'],loc='best')
savefig('1D_error.pdf')
show()
