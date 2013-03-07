# set9.py
# created 2/20/13 by stacy kim

import sys
from math import *
from numpy import *
from matplotlib.pyplot import *
from interpolate import *
from time import clock


data = array([[float(el) for el in line.split(' ')[1:] if el!='']
              for line in open('presupernova.dat','r').read().split('\n')[1:-1]])
r=data[:,1]
mass=data[:,0]
rho=data[:,3]

ROUT=1.0e9
G=6.67398e-8

"""
plot(r,rho)
xlabel('radius (cm)')
ylabel('density (g/cm^3)')
xscale('log')
yscale('log')
savefig('rho.pdf')
show()

# plot all over r
for i in range(len(data[0])-1):
    subplot(2,3,i+1)
    plot(data[:,1],data[:,i])
    xscale('log')
    if i != 4 and i!= 5: yscale('log')
savefig('all_columns_vs_r.pdf')
show()

# check if guesses to which columns are mass, radius, and rho are correct
dm=[4*pi/3*r[0]**3*rho[0]]+[4*pi/3*(r[i]**3-r[i-1]**3)*rho[i] for i in range(1,len(r))]
for mm in dm: 
    if mm<0: print mm
mass=[sum(dm[:i]) for i in range(len(r))]
for i in range(30): print r[i],mass[i]

plot(r,data[:,0],'r')
plot(r,mass)
xscale('log')
yscale('log')
savefig('mass_check.pdf')
show()
"""

# THE 1D POISSION EQUATION SOLVERS -----------------------------------------------

def rhs(r,x,rho):
    """Implements the 1D Poisson equation."""
    # x = [phi,z]
    dphi=x[1]
    dz=4*pi*G*rho - 2*x[1]/r

    return array([dphi,dz])


def directODE(r,rho,dr,rout):
    """
    Solves Poisson's equation via the direct ODE method for the given density
    profile r, rho on a grid that extends from 0 to rout.  If a constant density
    profile is to be assumed, a constant rho can be given (r will be ignored).

    The solution and the grid it was computed on (r_grid,phi) is returned.
    """
    # construct the grid and calculate Mtot
    r_grid=arange(0.,rout,dr)
    nr=len(r_grid)
    try:
        len(rho)
        rho_grid=pw_quadratic(r,rho,r_grid)

        for ir in range(len(rho)): 
            if r[ir] > ROUT: break
        M=4*pi/3*(rho[0]*r[0]**3 +
                  sum([rho[i]*(r[i]**3-r[i-1]**3) for i in range(1,ir)]))
    except TypeError:
        rho_grid=rho*ones(nr)
        M=4*pi/3*rout**3

    # solve Poisson's eqn
    t0=clock()
    x=[array([0,4*pi*G*rho_grid[0]])] # (phi,z) boundary conditions
    for i in range(1,nr):
        x.append(x[i-1]+dr*rhs(r_grid[i],x[i-1],rho_grid[i-1]))
    
    phi_offset=-G*M/rout-x[-1][0]
    phi=array([p+phi_offset for p,z in x])
    print 'dODE: took',clock()-t0,'sec'

    return r_grid,phi,phi_offset


def matrix_method(r,rho,dr,rout):
    """
    Solves Poisson's equation via the matrix method for the given density profile 
    r, rho on a grid that extends from dr/2 to rout+dr/2 (shifted to avoid the 
    singularity at r=0). If a constant density profile is to be assumed, a constant
    rho can be given (r will be ignored).

    The solution and the grid it was computed on (r_grid,phi) is returned.
    """
    # construct the grid and calculate Mtot
    r_grid=arange(0.,rout,dr)+0.5*dr
    nr=len(r_grid)
    try:
        len(rho)
        rho_grid=pw_quadratic(r,rho,r_grid)
        for ir in range(len(rho)): 
            if r[ir] > ROUT: break
        M=4*pi/3 * (r[0]**3*rho[0] +
                    sum([(r[i]**3-r[i-1]**3)*rho[i] for i in range(1,ir)]))
    except TypeError:
        rho_grid=rho*ones(nr)
        M=4*pi/3*ROUT**3

    # construct the system of equations to be solved
    t0=clock()
    J=zeros([nr,nr])
    J[0,0]=-1./dr**2-1/(r_grid[0]*dr)
    J[0,1]=1./dr**2+1./(r_grid[0]*dr)
    for j in range(1,nr):
        J[j,j-1] = 1./dr**2-1./(r_grid[j]*dr)
        J[j,j]   = -2./dr**2
        if j!=(nr-1): J[j,j+1] = 1./dr**2+1./(r_grid[j]*dr)

    b=4*pi*G*array(rho_grid)

    # solve Poisson's eqn
    ymatrix=linalg.solve(J,b)
    
    phi_offset=-G*M/rout-ymatrix[-1]
    phi=array([p+phi_offset for p in ymatrix])
    print 'mm: took',clock()-t0,'sec'

    return r_grid,phi,phi_offset


# VALIDATION --------------------------------------------------------------------
# Solve Poisson's eqn for homogeneous sphere
rho=1
drr=[1e6,1e7]

analytic = lambda r,rho: 2./3*pi*G*rho*(r**2 - 3*ROUT**2)
r_ana=arange(0,ROUT,1e6)
phi_ana=array([analytic(r,rho) for r in r_ana])


# via the direct ODE method
plot(r_ana,phi_ana)

err=[]
for i in range(len(drr)):
    r_grid,phi = directODE(0,rho,drr[i],ROUT)
    plot(r_grid,phi)
    err.append([phi[i]-analytic(r_grid[i],rho) for i in range(len(phi))])

xlabel('radius (cm)')
ylabel('gravitational potential')
legend(['analytic','numerical, dr=1e6','numerical, dr=1e7','numerical, dr=1e8'],loc='best')
savefig('validation_directODE.pdf')
show()

plot(arange(0,ROUT,drr[0]),err[0],
     arange(0,ROUT,drr[1]),err[1],
     arange(0,ROUT,drr[1]),array(err[1])/10.)
xlabel('radius (cm)')
ylabel('error')
legend(['dr=1e6','dr=1e7','error(dr=1e7)/10'],loc='best')
show()


# via the matrix method
plot(r_ana,phi_ana)

err=[]
for i in range(len(drr)):
    r_grid,phi = matrix_method(0,rho,drr[i],ROUT)
    plot(r_grid,phi)
    err.append([phi[i]-analytic(r_grid[i],rho) for i in range(len(phi))])

xlabel('radius (cm)')
ylabel('gravitational potential')
legend(['analytic','numerical, dr=1e6','numerical, dr=1e7','numerical, dr=1e8'],loc='best')
savefig('validation_matrix.pdf')
show()

plot(arange(0,ROUT,drr[0])+drr[0]/2,err[0],
     arange(0,ROUT,drr[1])+drr[1]/2,err[1],
     arange(0,ROUT,drr[1])+drr[0]/2,array(err[1])/10.)
xlabel('radius (cm)')
ylabel('error')
legend(['dr=1e6','dr=1e7','error(dr=1e7)/10'],loc='best')
show()


# THE PRESUPERNOVA MODEL --------------------------------------------------------
# Solve Poisson's eqn for given density profile 
r=data[:,1]
rho=data[:,3]
dr=1e6

r_directODE,phi_directODE,phi_off_dODE = directODE(r,rho,dr,ROUT) #direct ODE method
r_matrix,phi_matrix,phi_off_matrix = matrix_method(r,rho,dr,ROUT) # matrix method

plot(r_directODE,phi_directODE)
plot(r_matrix,phi_matrix)
xlabel('radius (cm)')
ylabel('gravitational potential')
xscale('log')
xlim([5e5,ROUT])
legend(['direct ODE method','matrix method'],loc='best')
savefig('comparison.pdf')
show()


# convergence tests
r_dODE=[r_directODE]
r_matrix=[r_matrix]
phi_dODE=[phi_directODE]
phi_matrix=[phi_matrix]
phi_off=[phi_off_dODE]

drr=[1e6,1e7,1e8]

plot(r_directODE,phi_directODE,'r')
plot(r_matrix,phi_matrix,'b')

for dr in drr[1:]:
    r_sol,phi,offset=directODE(r,rho,dr,ROUT)
    r_dODE.append(r_sol)
    phi_dODE.append(phi)
    phi_off.append(offset)
    plot(r_sol,phi,'r')

    r_sol,phi,offset=matrix_method(r,rho,dr,ROUT)
    r_matrix.append(r_sol)
    phi_matrix.append(phi)
    plot(r_sol,phi,'b')

show()

print phi_off

err_dODE=[[phi_dODE[i][j] - phi_dODE[0][where(r_dODE[0]==r_dODE[i][j])[0][0]]
           for j in range(len(r_dODE[i]))] for i in range(1,len(drr))]

err_dODE[0] = abs(array(err_dODE[0]) - phi_off[1] + phi_off[0])
err_dODE[1] = abs(array(err_dODE[1]) - phi_off[2] + phi_off[0])

plot(r_dODE[1],err_dODE[0],
     r_dODE[2],err_dODE[1],
     r_dODE[2],array(err_dODE[1])/10.)
xlabel('radius (cm)')
ylabel('error')
legend(['dr=1e7','dr=1e8','error(dr=1e8)/10'],loc='best')
savefig('dODE_preSN_convergence.pdf')
show()
