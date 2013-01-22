# set3.py
# created 1/18/13 by stacy kim
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from integrate import *


# EXERCISE 1 --------------------------------------------------------------------
# Integration via Newton-Cotes Formulae

# for the function sin(x)
trap=np.array([[math.pi/N, trapezoid(math.sin,0,math.pi,N)] for N in 2**np.arange(0,6)])
simp=np.array([[math.pi/N, simpson  (math.sin,0,math.pi,N)] for N in 2**np.arange(0,6)])

t_relerr=[abs(trap[i,1]-trap[i+1,1]) for i in range(len(trap)-1)]
s_relerr=[abs(simp[i,1]-simp[i+1,1]) for i in range(len(simp)-1)]
print 'sin(x)\nstep size         trapezoidal       simpson'
for i in range(len(trap)-2):
    print '{0:<17} {1:<17} {2:<17}'.format(trap[i,0],t_relerr[i]/t_relerr[i+1],
                                           s_relerr[i]/s_relerr[i+1])

plt.plot(trap[:,0],trap[:,1]-2.,'r',simp[:,0],simp[:,1]-2.,'b')
plt.legend(['trapezoidal','simpson'],loc='best')
plt.xlabel('Step Size')
plt.ylabel('Absolute Error')
plt.savefig('int_sin_convg.pdf')
plt.show()


# for the function x*sin(x)
xsin = lambda x: x*math.sin(x)
trap=np.array([[math.pi/N, trapezoid(xsin,0,math.pi,N)] for N in 2**np.arange(0,6)])
simp=np.array([[math.pi/N, simpson  (xsin,0,math.pi,N)] for N in 2**np.arange(0,6)])

t_relerr=[abs(trap[i,1]-trap[i+1,1]) for i in range(len(trap)-1)]
s_relerr=[abs(simp[i,1]-simp[i+1,1]) for i in range(len(simp)-1)]
print '\nx*sin(x)\nstep size         trapezoidal       simpson'
for i in range(len(trap)-2):
    print '{0:<17} {1:<17} {2:<17}'.format(trap[i,0],t_relerr[i]/t_relerr[i+1],
                                           s_relerr[i]/s_relerr[i+1])

plt.plot(trap[:,0],trap[:,1]-math.pi,'r',simp[:,0],simp[:,1]-math.pi,'b')
plt.legend(['trapezoidal','simpson'],loc='best')
plt.xlabel('Step Size')
plt.ylabel('Absolute Error')
plt.savefig('int_xsin_convg.pdf')
plt.show()


# EXERCISE 2 --------------------------------------------------------------------
# Gaussian Quadrature

EV  = 1.602177e-19 # in J
MEV = 1e6*EV       # in J
HBAR= 1.054571e-34 # in J s
C   = 2.99792e8    # in m/s

coeff = lambda kT:  math.pi*(kT/(math.pi*HBAR*C))**3
f = lambda x: x*x/(math.exp(x)+1)

# Gauss-Laguerre quadrature to compute electron number density
nodes=range(2,10)
integ=np.array([gauss_laguerre(f,0,float('inf'),n) for n in nodes])
ne   =coeff(20*MEV)*integ

print 'ne =',ne[0],ne[-1], ne[-2]/ne[-1]

relerr=[float(ne[i])/ne[-1]-1 for i in range(len(ne))]
plt.plot(nodes,relerr)
plt.xlim([1.5,nodes[-1]+0.5])
plt.xlabel('Number of nodes')
plt.ylabel('Relative Error')
plt.savefig('ne_convg.pdf')
plt.show()

# Gauss-Legendre quadrature to compute electron number density over E spectrum
ne_nodes=[]
E_range=np.arange(0,150,5.)
const=math.pi*(MEV/(math.pi*HBAR*C))**3
for n in nodes:
    integ_tot=0
    integ=[]
    for E in E_range:
        def f2(E):
            f1 = lambda x: x*x/(math.exp(x/E)+1) # where E in units of MEV
            return gauss_laguerre(f1,0,float('inf'),10)
        
        integ.append(gauss_legendre(f2,E,E+5,n))
        integ_tot+=integ[-1]

    if n==nodes[-1]:
        plt.plot(E_range,const*np.array(integ))
        plt.xlabel('Thermal energy (MeV)')
        plt.ylabel('Number density (gm^-1)')
        plt.savefig('ne_over_e_spectrum.pdf')
        plt.show()

    ne_nodes.append(integ_tot*const)

relerr=[float(ne_nodes[i])/ne_nodes[-1]-1 for i in range(len(ne))]
plt.plot(nodes, relerr)
plt.xlim([1.5,nodes[-1]+0.5])
plt.xlabel('Number of nodes')
plt.ylabel('Relative Error')
plt.savefig('ne_spectrum_convg.pdf')
plt.show()

print ne_nodes[-1], ne_nodes[-2]/ne_nodes[-1]

