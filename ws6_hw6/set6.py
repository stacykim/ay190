# set6.py
# created 1/29/13 by stacy kim

import numpy as np, time
import matplotlib.pyplot as plt
from gaussianelimination import *
from rungekutta import *


# EXERCISE 1 --------------------------------------------------------------------
print 'reading in LSEs with dimensions'
m,b=[],[]
for i in range(1,6):
    b.append([float(line[1:]) 
                  for line in open('LSE'+str(i)+'_bvec.dat').read().split('\n')[:-1]])
    m.append([[float(n) for n in line.split(' ')[1:] if n != ''] \
                  for line in open('LSE'+str(i)+'_m.dat').read().split('\n')[:-1]])

    print '  ',i, '\tm:',len(m[-1]), len(m[-1][0]),'\tb:',len(b[-1])


for i in range(len(m)):
    sy=[m[i][j]+[b[i][j]] for j in range(len(m[i]))]
    #t1=timeit.timeit(stmt='x1=myGauss(sy)',setup='from __main__ import myGauss, sy,x1',number=1)
    t1=time.clock()
    ans1=myGauss(sy)
    t1=time.clock()-t1

    t2=time.clock()
    ans2=np.linalg.solve(m[i],b[i])
    t2=time.clock()-t2
    print i,':',t1, t2

    if np.all([ans1[j]==ans2[j] for j in range(len(ans1))]):
        print 'Does not match for LSE',j+1,'!'


# EXERCISE 2 --------------------------------------------------------------------
# Stiff ODE Systems

# exact solution
t_exact =np.arange(0,4,0.01)
y1_exact=np.array([(99+math.exp(100*tt))/100 for tt in t_exact])
y2_exact=np.array([(1-math.exp(100*tt))/100 for tt in t_exact])

def plot_exact():
    plt.subplot(211)
    plt.plot(t_exact,y1_exact,'k',lw=2)
    plt.ylabel('Y1')
    plt.yscale('log')
    plt.xlim([0,4])

    plt.subplot(212)
    plt.plot(t_exact,-y2_exact,'k',lw=2)
    plt.xlabel('time')
    plt.ylabel('-Y2')
    plt.yscale('log')
    plt.xlim([0,4])


# RHS of our stiff ODE system
def rhs(x,t):
    # x = [y1, y2]
    dy1= x[0] - 99.0*x[1]
    dy2=-x[0] + 99.0*x[1]
    return np.array([dy1,dy2])


# Inputs
h=[0.1,0.01,1e-3,1e-4] # step size
x0=np.array([1.,0.])   # initial condition, x=[y1,y2]
t0=0.0                 # initial time
tf=4.0                 # time to integrate until
err=[1.0,1.0]          # desired accuracy for each element of x
ls=[':','-.','--','-']      # line styles (one for each each step size)


# explicit Euler
print 'explicit euler'
plot_exact()
for hh in h:
    x=[x0[0],x0[1]]
    y1,y2=[x[0]],[x[1]]
    t=np.arange(t0,tf,hh)
    for tt in t[1:]:
        x+=hh*rhs(x,tt)
        y1.append(x[0])
        y2.append(x[1])

    plt.subplot(211)
    plt.plot(t,y1,'r',ls=ls[h.index(hh)])
    plt.subplot(212)
    plt.plot(t,-np.array(y2),'r',ls=ls[h.index(hh)])

plt.savefig('stiff_explicit_euler.pdf')
plt.show()


# RK2 integrator
print 'rk2'
plot_exact()
for hh in h:
    fn='stiff_rk2_h{0}.dat'.format(hh)
    driver(rk2,x0,t0,tf,hh,rhs,err,fn)

    dat = np.array([[float(el) for el in line.split(' ') if el != ''] 
                    for line in open(fn).read().split('\n')[3:-1]])
    t,y1,y2=dat[:,0],dat[:,2],dat[:,3]

    plt.subplot(211)
    plt.plot(t,y1,'r',ls=ls[h.index(hh)])
    plt.subplot(212)
    plt.plot(t,-y2,'r',ls=ls[h.index(hh)])

plt.savefig('stiff_rk2.pdf')
plt.show()


# RK4 integrator
print 'rk4'
plot_exact()
for hh in h:
    fn='stiff_rk4_h{0}.dat'.format(hh)
    driver(rk4,x0,t0,tf,hh,rhs,err,fn)

    dat = np.array([[float(el) for el in line.split(' ') if el != ''] 
                    for line in open(fn).read().split('\n')[3:-1]])
    t,y1,y2=dat[:,0],dat[:,2],dat[:,3]

    plt.subplot(211)
    plt.plot(t,y1,'r',ls=ls[h.index(hh)])
    plt.subplot(212)
    plt.plot(t,-y2,'r',ls=ls[h.index(hh)])

plt.savefig('stiff_rk4.pdf')
plt.show()


# backward Euler
print 'bkwd euler'
plot_exact()
for hh in h:
    if hh==0.01: continue # will cause zero division error
    x=[x0[0],x0[1]]
    y1,y2=[x[0]],[x[1]]
    t=np.arange(t0,tf,hh)
    for tt in t[1:]:
        y1.append(-(y1[-1]-99*y1[-1]*hh-99*y2[-1]*hh)/(100*hh-1))
        y2.append(-(y2[-1]-y1[-1]*hh-y2[-1]*hh)/(100*hh-1))

    plt.subplot(211)
    plt.plot(t,y1,'r',ls=ls[h.index(hh)])
    plt.xlim([3.9,4])
    plt.ylim([1e167,1e182])
    plt.subplot(212)
    plt.xlim([3.9,4])
    plt.ylim([1e167,1e182])
    plt.plot(t,-np.array(y2),'r',ls=ls[h.index(hh)])

plt.savefig('stiff_bkwd_euler_closeup.pdf')
plt.show()
