############# equation of state

import sys,math
from pylab import *

def write_state(hyd,fn):
    f=open(fn,'w')
    f.write('  i          x            density        pressure        velocity\n\n')
    for j in range(hyd.n):
        f.write('{4:4d} {0:15E} {1:15E} {2:15E} {3:15E}\n'.format(hyd.x[j],hyd.rho[j],hyd.press[j],hyd.vel[j],j))
    f.close()


# ================================================================================
# For stellar core collapse setup

K1 = 1.2435e15 * (0.5e0**(4.0/3.0))
gamma1 = 1.28
gamma2 = 2.5
gammath = 1.5
rhonuc = 2.0e14

E1 = K1/(gamma1-1.e0)
E2 = (gamma1 - 1.e0)/(gamma2-1.e0)*E1*rhonuc**(gamma1-gamma2)
K2 = (gamma2 - 1.e0)*E2
E3 = (gamma2 - gamma1)/(gamma2-1.e0)*E1*rhonuc**(gamma1-1.e0)


def hybrid_eos(rho,eps,hyd):
    press = zeros(len(rho))
    cs2 = zeros(len(rho))

    for i in range(len(rho)):
        if(rho[i] < rhonuc):
            Kx=K1
            Ex=E1
            Gx=gamma1
            Ex3=0.e0
        else:
            Kx=K2
            Ex=E2
            Gx=gamma2
            Ex3=E3

        up=Ex*rho[i]**Gx+Ex3*rho[i]
        eth = eps[i]*rho[i] - up
        pth = (gammath - 1)*eth
        dpth_drho=(gammath - 1)*(eps[i]-Ex*Gx*rho[i]**(Gx-1.e0)-Ex3)
        dpth_deps=(gammath - 1)*rho[i]
        
        pco = Kx*rho[i]**Gx
        dpco_drho=Gx*pco/(rho[i]+1.0e-20)
        dpco_deps=0.e0
        
        press[i] =pco+pth
        dpde = dpco_deps+dpth_deps
        dp_drho = dpco_drho+dpth_drho
        cs2[i] = dp_drho+dpde*press[i]/(rho[i]+1.0e-20)**2

        if cs2[i] < 0 and i >= hyd.g and i <= (hyd.n - hyd.g):
            print i,dp_drho,dpde,press[i],rho[i]
            write_state(hyd,'error.dat')
            #sys.exit()
            
    return (press,cs2)


def hybrid_main():

    rho=2.0e14
    eps=2.0e20

    (press,cs2) = hybrid_eos(rho,eps)
    
    print "%15.6E %15.6E" % (press,sqrt(cs2))


# ================================================================================
# For shocktube setup

def eos_press(rho,eps,gamma):
    press = (gamma - 1.0) * rho * eps
    return press

def eos_cs2(rho,eps,gamma):
    prs = (gamma - 1.0) * rho *eps
    dpde = (gamma - 1.0) * rho
    dpdrho = (gamma - 1.0) * eps
    cs2 = dpdrho + dpde * prs/(rho+1.0e-30)**2
    if (cs2 < 0).any(): print 'rho =',rho,'\neps =',eps
    return cs2
