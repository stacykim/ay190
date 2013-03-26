#!/usr/bin/env python
import sys,math
from pylab import *
from eos import *

# basic parameters
K1 = 1.2435e15 * (0.5e0**(4.0/3.0))
gamma1 = 1.28
gamma2 = 2.5
gammath = 1.5
rhonuc = 2.0e14

E1 = K1/(gamma1-1.e0)
E2 = (gamma1 - 1.e0)/(gamma2-1.e0)*E1*rhonuc**(gamma1-gamma2)
K2 = (gamma2 - 1.e0)*E2
E3 = (gamma2 - gamma1)/(gamma2-1.e0)*E1*rhonuc**(gamma1-1.e0)

G = 6.672e-8 # gravitational constant
rhomin=1e6

cfl = 0.5
dt = 1.0e-5
dtp = dt
reconstruction_type = 'pc' # pc, minmod, mc
nzones = 2000
tend = 1.0

##################################### class definition
class mydata:
    def __init__(self,nzones):
        self.x      = zeros(nzones) # cell centers
        self.xi     = zeros(nzones) # cell LEFT interfaces
        self.rho    = zeros(nzones)
        self.rhop   = zeros(nzones)
        self.rhom   = zeros(nzones)
        self.vel    = zeros(nzones)
        self.velp   = zeros(nzones)
        self.velm   = zeros(nzones)
        self.eps    = zeros(nzones)
        self.epsp   = zeros(nzones)
        self.epsm   = zeros(nzones)
        self.press  = zeros(nzones)
        self.pressp = zeros(nzones)
        self.pressm = zeros(nzones)
        self.q      = zeros((3,nzones)) # conserved quantities
        self.qp     = zeros((3,nzones))  
        self.qm     = zeros((3,nzones))  
        self.n      = nzones
        self.g      = 3 # ghost cells at each of inner/outer edges of domain


    def setup_grid(self,xmin,xmax):
        dx = (xmax - xmin) / (self.n - self.g*2 - 1)
        xmin = xmin - self.g*dx
        xmax = xmax + self.g*dx
        for i in range(self.n):
            self.x[i] = xmin + (i)*dx        # cell centers
            self.xi[i] = self.x[i] - 0.5*dx  # cell LEFT interfaces            


##################################### initial conditions
def setup_star(hyd):
    poly = loadtxt('poly.dat')
    prad = poly[:,0]
    prho = poly[:,1]
    nn = len(prho)

    for i in range(hyd.n):
        hyd.rho[i] = max(rhomin,linterp(hyd.x[i],nn,prho,prad))
        hyd.press[i] = K1 * hyd.rho[i]**gamma1
        hyd.eps[i] = hyd.press[i] / (gamma1-1.0) / hyd.rho[i]
        if(hyd.rho[i] <= rhomin):
            hyd.rho[i] = hyd.rho[i] / 5.0
    
    return hyd


def linterp(xx,n,f,x):
    i=0
    while(i<n and x[i] < xx):  i=i+1

    if (i==n):    ff = rhomin
    elif (i==0):  ff = (f[1]-f[0])/(x[1]-x[0]) * (xx - x[0]) + f[0]
    else:         ff = (f[i]-f[i-1])/(x[i]-x[i-1]) * (xx-x[i-1]) + f[i-1]
        
    return ff


##################################### some basic functions
def prim2con(rho,vel,eps):
    q = zeros((3,len(rho)))
    q[0,:] = rho[:]
    q[1,:] = rho[:] * vel[:]
    q[2,:] = rho[:] * eps[:] + 0.5 * rho[:] * vel[:]**2

    return q

def con2prim(q):
    global hyd

    rho = q[0,:]
    vel = q[1,:] / rho
    eps = q[2,:] / rho - 0.5*vel**2
    press,cs2 = hybrid_eos(rho,eps,hyd)

    return (rho,eps,press,vel)

############# minmod function 
def minmod(a,b):
    if(a*b < 0):
        mm = 0.0
    elif(abs(a)<abs(b)):
        mm=a
    else:
        mm=b
    return mm

############# minmod function 
def tvd_minmod_reconstruct(n,g,f,x,xi):
    fp = zeros(n)
    fm = zeros(n)
    for i in range(g-1,n-g+1):
        dx_up = x[i] - x[i-1]
        dx_down = x[i+1] - x[i]
        dx_m = x[i] - xi[i]
        dx_p = xi[i+1] - x[i]
        df_up = (f[i]-f[i-1]) / dx_up
        df_down = (f[i+1]-f[i]) / dx_down
        delta = minmod(df_up,df_down)
        fp[i] = f[i] + delta*dx_p
        fm[i] = f[i] - delta*dx_m

    return (fp,fm)

############# signum functions
def signum(x,y):
    if(y >= 0):
        return abs(x)
    else:
        return -abs(x)

############# mc reconstruction
def tvd_mc_reconstruct(n,g,f,x,xi):
    fp = zeros(n)
    fm = zeros(n)
    for i in range(g-1,n-g+1):
        dx_up = x[i] - x[i-1]
        dx_down = x[i+1] - x[i]
        dx_m = x[i] - xi[i]
        dx_p = xi[i+1] - x[i]
        df_up = (f[i]-f[i-1]) / dx_up
        df_down = (f[i+1]-f[i]) / dx_down
        if(df_up*df_down < 0):
            delta = 0.0
        else:
            delta = signum(min(2.0*abs(df_up),2.0*abs(df_down),\
                               0.5*(abs(df_up)+abs(df_down))),\
                               df_up + df_down)

        fp[i] = f[i] + delta*dx_p
        fm[i] = f[i] - delta*dx_m

    return (fp,fm)

############# reconstruction top level function
def reconstruct(hyd,type):
    if(type=='pc'):
        # piecewise constant reconstruction 
        for i in range(hyd.g-1,hyd.n-hyd.g+1):
            hyd.rhop[i] = hyd.rho[i]
            hyd.rhom[i] = hyd.rho[i]
            hyd.epsp[i] = hyd.eps[i]
            hyd.epsm[i] = hyd.eps[i]
            hyd.velp[i] = hyd.vel[i]
            hyd.velm[i] = hyd.vel[i]


    elif(type=='minmod'):
        (hyd.rhop,hyd.rhom) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.rho,hyd.x,hyd.xi)
        (hyd.epsp,hyd.epsm) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.eps,hyd.x,hyd.xi)
        (hyd.velp,hyd.velm) = tvd_minmod_reconstruct(hyd.n,hyd.g,hyd.vel,hyd.x,hyd.xi)

    elif(type=='mc'):
        (hyd.rhop,hyd.rhom) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.rho,hyd.x,hyd.xi)
        (hyd.epsp,hyd.epsm) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.eps,hyd.x,hyd.xi)
        (hyd.velp,hyd.velm) = tvd_mc_reconstruct(hyd.n,hyd.g,hyd.vel,hyd.x,hyd.xi)
                
    else:
        print "reconstruction type not known; abort!"
        sys.exit()


    hyd.pressp,cs2p = hybrid_eos(hyd.rhop,hyd.epsp,hyd)
    hyd.pressm,cs2m = hybrid_eos(hyd.rhom,hyd.epsm,hyd)

    hyd.qp = prim2con(hyd.rhop,hyd.velp,hyd.epsp)
    hyd.qm = prim2con(hyd.rhom,hyd.velm,hyd.epsm)

    return hyd


############# time step calculation
def calc_dt(hyd,dtp):
    global i
    press,cs2 = hybrid_eos(hyd.rho,hyd.eps,hyd)
    cs = sqrt(cs2)
    dtnew = [(hyd.x[j+1]-hyd.x[j]) / max(abs(hyd.vel[j]+cs[j]), abs(hyd.vel[j]-cs[j]))
             for j in range(hyd.g,hyd.n-hyd.g)]
    dtnew = min(dtnew)
    dtnew = min(cfl*dtnew,1.05*dtp)
    return dtnew

############# HLLE solver
def hlle(hyd):
    # compute eigenvalues
    evl  = zeros((3,hyd.n))
    evr  = zeros((3,hyd.n))
    smin = zeros(hyd.n)
    smax = zeros(hyd.n)

    pressp,cs2p = hybrid_eos(hyd.rhop,hyd.epsm,hyd)
    pressm,cs2m = hybrid_eos(hyd.rhop,hyd.epsm,hyd)
    csp  = sqrt(cs2p)
    csm  = sqrt(cs2m)

    for i in range(1,hyd.n-2):
        evl[0,i] = hyd.velp[i] - csp[i]
        evl[1,i] = hyd.velp[i]
        evl[2,i] = hyd.velp[i] + csp[i]

        evr[0,i] = hyd.velm[i+1] - csm[i+1]
        evr[1,i] = hyd.velm[i+1]
        evr[2,i] = hyd.velm[i+1] + csm[i+1]
        
        smin[i] = min(concatenate((evl[:,i],evr[:,i])))
        smax[i] = max(concatenate((evl[:,i],evr[:,i])))

    # set up flux left L and right R of the interface
    # at i+1/2
    fluxl = zeros((3,hyd.n))
    fluxr = zeros((3,hyd.n))
    
    for i in range(1,hyd.n-2): # skip only 1 boundary cell
        fluxl[0,i]=hyd.rhop[i] * hyd.velp[i]
        fluxl[1,i]=hyd.rhop[i] * hyd.velp[i]**2 + hyd.pressp[i]
        fluxl[2,i]=(hyd.rhop[i] * hyd.epsp[i] +
                    0.5*hyd.rhop[i] * hyd.velp[i]**2 +
                    hyd.pressp[i]) * hyd.velp[i]

        fluxr[0,i]=hyd.rhom[i+1] * hyd.velm[i+1]
        fluxr[1,i]=hyd.rhom[i+1] * hyd.velm[i+1]**2 + hyd.pressm[i+1]
        fluxr[2,i]=(hyd.rhom[i+1] * hyd.epsm[i+1] +
                    0.5*hyd.rhom[i+1] * hyd.velm[i+1]**2 +
                    hyd.pressm[i+1]) * hyd.velm[i+1]

    # solve the Riemann problem for the i+1/2 interface
    ds = smax - smin
    flux = zeros((3,hyd.n))
    for i in range(hyd.g-1,hyd.n-hyd.g+1):
        for j in range(3):
            flux[j,i] = (smax[i]*fluxl[j,i] - smin[i]*fluxr[j,i] + 
                         smin[i]*smax[i] * (hyd.qm[j,i+1]-hyd.qp[j,i])) / (smax[i]-smin[i])

    # flux differences
    fluxdiff = zeros((3,hyd.n))
    for i in range(hyd.g,hyd.n-hyd.g):
        dr = hyd.xi[i+1]-hyd.xi[i]
        for j in range(3):
            fluxdiff[j,i] = (hyd.xi[i+1]**2*flux[j,i] - hyd.xi[i]**2*flux[j,i-1]) \
                / (dr*hyd.xi[i+1]**2)

    return fluxdiff

############# mass calculation (required for RHS calcluation)
def calc_mass(hyd):
    mass = zeros(hyd.n)
    m1 = zeros(hyd.n)
    m1[hyd.g] = 4.0/3.0*pi * hyd.rho[0] * (hyd.xi[hyd.g+1]**3)
    mass[hyd.g] = 4.0/3.0*pi * hyd.rho[0] * (hyd.x[hyd.g])**3
    for i in range(hyd.g,hyd.n-1):
        m1[i] = 4.0/3.0*pi * hyd.rho[i] * (hyd.xi[i+1]**3 - hyd.xi[i]**3)
        dmi = 4.0/3.0*pi * (hyd.rho[i-1]*(hyd.xi[i]**3 - hyd.x[i-1]**3) +
                            hyd.rho[i]*(hyd.x[i]**3 - hyd.xi[i]**3))
        mass[i] = mass[i-1] + dmi
        
    return (mass,m1)


############# RHS calculation
def calc_rhs(hyd):
    # reconstruction and prim2con
    hyd = reconstruct(hyd,reconstruction_type)
    # compute flux differences
    fluxdiff = hlle(hyd)

    # compute non-conserved terms
    rhs = zeros((3,hyd.n))
    mass,m1 = calc_mass(hyd)
    for i in range(hyd.g,hyd.n-hyd.g):
        dr = hyd.xi[i+1]-hyd.xi[i]
        #dr = (hyd.xi[i+1]**3 - hyd.xi[i]**3) / (3*hyd.x[i]**2)
        rhs[0,i] = -fluxdiff[0,i]
        rhs[1,i] = (-fluxdiff[1,i]
                     - hyd.rho[i+1]*G*(m1[i+1]-m1[i])/(hyd.xi[i+1]**2 * dr)
                     - (hyd.press[i+1]-hyd.press[i])/dr)
        rhs[2,i] = (-fluxdiff[2,i]
                     - hyd.vel[i]*hyd.rho[i]*G*(m1[i+1]-m1[i])/(hyd.xi[i+1]**2 *dr))

    return rhs

############# boundary conditions
def apply_bcs_spherical(hyd):
    # hardcoded for 3 inner boundary points
    hyd.rho[0] = hyd.rho[5]
    hyd.rho[1] = hyd.rho[4]
    hyd.rho[2] = hyd.rho[3]
    hyd.vel[0] = -hyd.vel[5]
    hyd.vel[1] = -hyd.vel[4]
    hyd.vel[2] = -hyd.vel[3]
    hyd.eps[0] = hyd.eps[5]
    hyd.eps[1] = hyd.eps[4]
    hyd.eps[2] = hyd.eps[3]
    hyd.press[0] = hyd.press[5]
    hyd.press[1] = hyd.press[4]
    hyd.press[2] = hyd.press[3]

    hyd.rho[hyd.n-hyd.g:hyd.n-1] = hyd.rho[hyd.n-hyd.g-1]
    hyd.vel[hyd.n-hyd.g:hyd.n-1] = hyd.vel[hyd.n-hyd.g-1]
    hyd.eps[hyd.n-hyd.g:hyd.n-1] = hyd.eps[hyd.n-hyd.g-1]
    hyd.press[hyd.n-hyd.g:hyd.n-1] = hyd.press[hyd.n-hyd.g-1]

    return hyd

def write_state(hyd,fn):
    f=open(fn,'w')
    f.write('  i          x            density        pressure        velocity\n\n')
    for j in range(hyd.n):
        f.write('{4:4d} {0:15E} {1:15E} {2:15E} {3:15E}\n'.format(hyd.x[j],hyd.rho[j],hyd.press[j],hyd.vel[j],j))
    f.close()


########################################################################
# Main program
########################################################################


hyd = mydata(nzones)     # initialize
hyd.setup_grid(0.0,2e9)  # set up grid
setup_star(hyd)          # set up initial data
dt = calc_dt(hyd,dt)     # get initial timestep

velold=0
write_state(hyd,'collapse_0ms.dat')
frho = open('central_density.dat','w')

# initial prim2con
hyd.q = prim2con(hyd.rho,hyd.vel,hyd.eps)

t = 0.0
tcb = 9e99 # time of core bounce
i = 0

passed5ms  = 0
passed_core_bounce = 0
core_bounce10ms = 0

# display stuff
figure()
ion()
show()
max_rho0 = max(hyd.rho)
max_press0 = max(hyd.press)
print 'initial maxima:  rho =',max_rho0,'pressure =',max_press0

# main integration loop
while t < tend:

    frho.write('{0:15.6E} {1:15.6E}\n'.format(t,hyd.rho[hyd.g]))

    if i % 10 == 0 or (t > 0.005 and not passed5ms) \
            or (hyd.vel[hyd.g] > 0 and i >= 10 and not passed_core_bounce) \
            or (t > (tcb+0.01) and not core_bounce10ms) \
            or (t > 1.99e-2):

        clf()
        plot(hyd.x,hyd.rho/max_rho0,'b')
        plot(hyd.x,hyd.press/max_press0,'r')
        plot(hyd.x,hyd.vel/1e9,'g')
        xlim([1e6,2e9])
        xlabel('radius (cm)')
        legend(['density','pressure','velocity'],loc='best')        
        xscale('log')

        msg = ''

        if (t > 0.005 and not passed5ms):
            savefig('collapse_5ms.pdf')
            write_state(hyd,'collapse_5ms.dat')
            msg = '5ms'
            passed5ms = 1

        if (hyd.vel[hyd.g] > 0 and i >= 10 and not passed_core_bounce):
            savefig('collapse_corebounce.pdf')
            write_state(hyd,'collapse_corebounce.dat')
            tcb = t
            msg = 'cb'
            passed_core_bounce = 1

        if (t >= (tcb+0.01) and not core_bounce10ms):
            savefig('collapse_corebounce10ms.pdf')
            write_state(hyd,'collapse_corebounce10ms.dat')
            msg = 'cb10ms'
            core_bounce10ms = 1

        if (t > 1.99e-2):
            ylim([-2,26])            
            savefig('collapse_last.pdf')
            write_state(hyd,'collapse_last.dat')
            msg = 'last'

        print "%5d %15.6E %15.6E %15.6E %s" % (i,t,dt,hyd.rho[hyd.g],msg)

        #ylim([-1.5,100])        
        draw()

    # calculate new timestep
    dt = calc_dt(hyd,dt)

    # save old state
    hydold = hyd
    velold = hyd.vel[hyd.g]
    qold = hyd.q

    # calc rhs
    k1 = calc_rhs(hyd)
    # calculate intermediate step
    hyd.q = qold + 1.0/2.0 * dt * k1
    # con2prim
    (hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q)
    # boundaries
    hyd = apply_bcs_spherical(hyd)

    #calc rhs
    k2 = calc_rhs(hyd)
    #apply update
    hyd.q = qold + dt * 0.5 * (k1 + k2)
    # con2prim
    (hyd.rho,hyd.eps,hyd.press,hyd.vel) = con2prim(hyd.q)
    # apply bcs
    hyd = apply_bcs_spherical(hyd)

    # update time
    t += dt
    i += 1

# display final result
ioff()
legend(['density','pressure','velocity'],loc='best')
savefig('hlle_'+reconstruction_type+'.pdf')
show()

