# set2.py
# created 1/14/12 by stacy kim

import interp
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Cephieid Lighcurve data
xdat=[0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0]
ydat=[0.302,0.185,0.106,0.093,0.24,0.579,0.561,0.468,0.302]


# IN-CLASS EXERCISE 1 -----------------------------------------------------------
# Finite difference approximation and convergence
h=1.
x1=np.arange(-2,6,h)
f1=x1**3-5*x1*x1+x1

# calculate first derivative
fp1=3*x1*x1-10*x1+1                                          # exact
fp_fd1=[(f1[i+1]-f1[i])/h for i in range(len(x1)-1)]         # forward diff, h
fp_cd1=[(f1[i+1]-f1[i-1])/(2*h) for i in range(1,len(x1)-1)] # central diff, h

# redo with reduced step size
h=.5
x2=np.arange(-2,6,h)
f2=x2**3-5*x2*x2+x2

fp2=3*x2*x2-10*x2+1                                          # exact
fp_fd2=[(f2[i+1]-f2[i])/h for i in range(len(x2)-1)]         # forward diff, h/2
fp_cd2=[(f2[i+1]-f2[i-1])/(2*h) for i in range(1,len(x2)-1)] # central diff, h/2

# plot approximations
x=np.arange(-2,6,0.05)
fp=3*x*x-10*x+1

plt.plot(x,fp,c='k')
plt.plot(x1[:-1],fp_fd1,c='r',ls=':')
plt.plot(x2[:-1],fp_fd2,c='r',ls='--')
plt.savefig('forward_diff.pdf')
plt.show()

plt.plot(x,fp,c='k')
plt.plot(x1[1:-1],fp_cd1,c='b',ls=':')
plt.plot(x2[1:-1],fp_cd2,c='b',ls='--')
plt.savefig('central_diff.pdf')
plt.show()

# convergence plots
plt.plot(x1[:-1],abs(fp_fd1-fp1[:-1]),c='r',ls=':')
plt.plot(x2[:-1],abs(fp_fd2-fp2[:-1]),c='r',ls='--')
plt.savefig('forward_diff_converg.pdf')
plt.show()

plt.plot(x1[1:-1],abs(fp_cd1-fp1[1:-1]),c='b',ls=':')
plt.plot(x2[1:-1],abs(fp_cd2-fp2[1:-1]),c='b',ls='--')
plt.savefig('central_diff_converg.pdf')
plt.show()


# IN-CLASS EXERCISE 2a ----------------------------------------------------------
# Lagrange interpolation of Cepheid lightcurve
x=np.arange(0,1,0.001)
y_lgr=interp.lagrange(xdat,ydat,x)

# plot interpolated function and original data
plt.plot(xdat,ydat,'ko',mec='k')
plt.plot(x,y_lgr,c='r')
plt.xlim([-0.1,1.1])
plt.xlabel('Time')
plt.ylabel('Apparent Magnitude')
plt.savefig('cepheid_lagrange.pdf')
plt.show()

# IN-CLASS EXERCISE 2b ----------------------------------------------------------
y_pwl=interp.pw_linear(xdat,ydat,x) # piecewise linear interp
y_pwq=interp.pw_quadratic(xdat,ydat,x) # piecewise quadratic

plt.plot(xdat,ydat,'ko',mec='k')
plt.plot(x,y_pwl,c='r')
plt.plot(x,y_pwq,c='b')
plt.legend(['data','piecewise linear','piecewise quadratic'],loc='best')
plt.xlim([-0.1,1.1])
plt.xlabel('Time')
plt.ylabel('Apparent Magnitude')
plt.savefig('cepheid_pwl_pwq.pdf')
plt.show()


# HOMEWORK EXERCISE 3 -----------------------------------------------------------
# piecewise cubic Hermite
x_pwch=np.arange(0.2,0.8,0.001)
ypdat=[(ydat[i+1]-ydat[i-1])/(xdat[i+1]-xdat[i-1]) for i in range(1,len(xdat)-1)]
y_pwch=interp.cubic_hermite(xdat[1:-1],ydat[1:-1],ypdat,x_pwch) 

# cubic spline
f_cs=interp1d(xdat,ydat,kind='cubic')

plt.plot(xdat,ydat,'ko',mec='k')
plt.plot(x_pwch,y_pwch,c='r')
plt.plot(x,f_cs(x),c='b')
plt.legend(['data','piecewise cubic Hermite','cubic spline'],loc='best')
plt.xlim([-0.1,1.1])
plt.xlabel('Time')
plt.ylabel('Apparent Magnitude')
plt.savefig('cepheid_pwch_cs.pdf')
plt.show()


# plot all results
plt.plot(xdat,ydat,'ko',mec='k')
plt.plot(x,y_lgr, x,y_pwl, x,y_pwq, x_pwch,y_pwch, x,f_cs(x))
plt.legend(['data','lagrange','piecewise linear','piecewise quadratic','cubic hermite','cublic spline'],loc='best')
plt.xlim([-0.1,1.1])
plt.ylim([0,1])
plt.xlabel('Time')
plt.ylabel('Apparent Magnitude')
plt.savefig('cepheid_all.pdf')
plt.show()
