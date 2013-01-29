# set5.py
# created 1/23/13 by stacy kim

import numpy as np, math, sys
import matplotlib.pyplot as plt
from findroot import *

# EXERCISE 1 --------------------------------------------------------------------
# Root Finding: Eccentricity Anomality

# orbital parameters of Earth
e=0.0167
a=1.496e6 # km
b=a*math.sqrt(1-e*e)
P=365.25635
omega=2*math.pi/P # angular velocity
t=91.0

f  = lambda E: E-omega*t+e*math.sin(E)
df = lambda E: 1+e*math.cos(E)

# plot the function at t=91.0 d
EE=np.arange(0,2*math.pi,2*math.pi/1e3)
plt.plot(EE,[f(ee) for ee in EE],'r')
t=181.0
plt.plot(EE,[f(ee) for ee in EE],'b')
t=273.0
plt.plot(EE,[f(ee) for ee in EE],'g')
plt.legend(['t = 91 days','t = 181 days','t = 273 days'],loc='best')
plt.plot([0,2*math.pi],[0,0],'k')
plt.xlim([EE[0],EE[-1]])
plt.xlabel('eccentric anomaly')
plt.ylabel('f(E)')
plt.savefig('f(91d).pdf')
plt.show()

# solve for position of Earth at 3 timepoints
tt=[91.0, 182.0, 273.0]
E0,err0,cnt0=np.zeros([3,len(tt),3])
for i in range(len(tt)):
    t=tt[i]
    E0[i,0],err0[i,0],cnt0[i,0]=newton_raphson(i*math.pi/2,f,df,1e-10)
    E0[i,1],err0[i,1],cnt0[i,1]=secant(0.95*(i+1)*math.pi/2,1.05*(i+1)*math.pi/2,f,1e-10)
    E0[i,2],cnt0[i,2]=bisection(0.95*(i+1)*math.pi/2,1.05*(i+1)*math.pi/2,f,1e-10)

x=np.array([[a*math.cos(ee) for ee in time] for time in E0])
y=np.array([[b*math.sin(ee) for ee in time] for time in E0])

# solve for position of Earth at many timepoints
full_orbit=[]
for i in np.arange(0,P,P/1e3):
    t=i
    full_orbit.append(secant(-1,7,f,1e-10)[0])
    
full_x=[a*math.cos(ee) for ee in full_orbit]
full_y=[b*math.sin(ee) for ee in full_orbit]

# plot resulting positions
plt.plot(x[:,0],y[:,0],'ro',x[:,1],y[:,1],'bo',x[:,2],y[:,2],'go')
plt.legend(['newton-rhapson','secant','bisection'],loc='best')
plt.plot(full_x,full_y,'k')
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.savefig('earth_pos.pdf')
plt.show()

# plot number of iterations required for solution
plt.plot(tt,cnt0[:,0],'r',tt,cnt0[:,1],'b',tt,cnt0[:,2],'g')
plt.legend(['newton-rhapson','secant','bisection'],loc='best')
plt.xlabel('time (days)')
plt.ylabel('number of iterations')
plt.savefig('niter_orig.pdf')
plt.show()



# HIGHLY ECCENTRIC EARTH ORBIT
e=0.99999
b=a*math.sqrt(1-e*e)

# plot the function at all t
EE=np.arange(0,2*math.pi,2*math.pi/1e3)
t=91.0
plt.plot(EE,[f(ee) for ee in EE],'r')
t=182.0
plt.plot(EE,[f(ee) for ee in EE],'b')
t=273.0
plt.plot(EE,[f(ee) for ee in EE],'g')
plt.legend(['t = 91 days','t = 182 days','t = 273 days'],loc='best')
plt.plot([0,2*math.pi],[0,0],'k')
plt.xlim([EE[0],EE[-1]])
plt.xlabel('eccentric anomaly')
plt.ylabel('f(E)')
plt.savefig('f(91d)_eccentric.pdf')
plt.show()

# solve for position of Earth at 3 timepoints
tt=[91.0, 182.0, 273.0]
E1,err1,cnt1=np.zeros([3,len(tt),3])
bounds=np.array([[0.0,1.0],[2.0,3.0],[5.0,6.0]])

for i in range(len(tt)):
    t=tt[i]
    E1[i,0],err1[i,0],cnt1[i,0]=newton_raphson(i*math.pi/2,f,df,1e-10)
    E1[i,1],err1[i,1],cnt1[i,1]=secant(bounds[i,0],bounds[i,1],f,1e-10)
    E1[i,2],cnt1[i,2]=bisection(bounds[i,0],bounds[i,1],f,1e-10)

x=np.array([[a*math.cos(ee) for ee in time] for time in E1])
y=np.array([[b*math.sin(ee) for ee in time] for time in E1])

# solve for position of Earth at many timepoints
full_orbit=[]
for i in np.arange(0,P,P/1e3):
    t=i
    full_orbit.append(secant(-1,7,f,1e-10)[0])
full_x=[a*math.cos(ee) for ee in full_orbit]
full_y=[b*math.sin(ee) for ee in full_orbit]

# plot resulting positions
plt.plot(x[:,0],y[:,0],'ro',x[:,1],y[:,1],'bo',x[:,2],y[:,2],'go')
plt.legend(['newton-rhapson','secant','bisection'],loc='best')
plt.plot(full_x,full_y,'k')
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.savefig('earth_pos_eccentric.pdf')
plt.show()

# plot number of iterations required for solution
plt.plot(tt,cnt1[:,0],'r',tt,cnt1[:,1],'b',tt,cnt1[:,2],'g')
plt.plot(tt,cnt0[:,0],'r--',tt,cnt0[:,1],'b--',tt,cnt0[:,2],'g--')
plt.legend(['newton-rhapson','secant','bisection'],loc='best')
plt.xlabel('time (days)')
plt.ylabel('number of iterations')
plt.savefig('niter_orig_eccentric.pdf')
plt.show()

# print number of iterations required for solution (realisitc & eccentric orbits)
print 'iterations required'
print 't= 91d\t', cnt0[:,0],'\n      \t',cnt1[:,0]
print 't=181d\t', cnt0[:,1],'\n      \t',cnt1[:,1]
print 't=273d\t', cnt0[:,2],'\n      \t',cnt1[:,2]

# EXERCISE 2 --------------------------------------------------------------------
# Root Finding: Polynomials with Multiple Roots
coeff=[0,0,0,-1,5,3] # f = sum([coeff[i] * x**i for i=0..len(coeff)])
f = lambda x: sum([coeff[i] * x**i for i in range(len(coeff))])

# plot the function
x=np.arange(-2,1,0.01)
y=[f(xx) for xx in x]
plt.plot(x,y,[x[0],x[-1]],[0,0],'k',[0,0],[-1,6],'k')
plt.xlim([x[0],x[-1]])
plt.ylim([-1,6])
plt.xlabel('x')
plt.ylabel('f(x)')
plt.savefig('polynomial_plot.pdf')
plt.show()

# calclate all roots of the polynomial
print 'roots =',find_all_roots(coeff)

