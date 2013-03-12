# set15.py
# created 3/12/13 by stacy kim

from pylab import *
from numpy import *
import sys


# MEASURING PI WITH AN MC EXPERIMENT ---------------------------------------------
niter=1000
nN=5
results=[]
ion()
for i in range(nN):
    N=10**i
    tot=0
    xx=[]
    yy=[]
    subplot(5,5,i+1)
    xlim([0,niter])
    ylim([1.5,4.1])
    for j in range(niter):
        random.seed(j)
        pairs = random.rand(N,2)*2-1 # x,y both \el [-1,1]
        n=sum([1 if (x*x+y*y)**0.5 < 1.0 else 0 for x,y in pairs])
        tot+=n
        if j % 50 == 0:
            xx+=[j]
            yy+=[4.*n/N]
            plot(xx,yy,'r+')
            draw()
    est_pi=4.*tot/N/niter
    print 'average estimate after',niter,'trials of',N,'samplings:',est_pi
    results.append(est_pi)

ioff()
savefig('sample_pi.pdf')
show()

clf()
x=pow(10,arange(0,5,0.1))
plot(x,[1/sqrt(val) for val in x],'k')
plot([10**i for i in range(nN)],abs(array(results)-pi))
xlabel('samples')
ylabel('error')
xscale('log')
yscale('log')
savefig('pi_mc_convg.pdf')
show()

"""
# RANDOM WALK (DIFFUSION) IN 1D --------------------------------------------------

step_size=0.01
all_steps=[]
niter=10000
ion()
xlim([0,niter])
xlabel('iteration')
ylabel('number of steps')
for i in range(niter):
    random.seed(i)
    steps=0
    path=0
    while abs(path) < 1:
        path+=step_size*(-1 if random.rand() < 0.5 else 1)
        steps+=1
    all_steps+=[steps]
    if i % 100 == 0 and i >= 100:
        plot(range(0,i+100,100),[all_steps[j] for j in range(0,i+100,100)],'r+')
        draw()

ioff()
savefig('walk_1D.pdf')
show()
print sum(all_steps)/float(niter)


# RANDOM WALK (DIFFUSION) IN 2D --------------------------------------------------

stepsize=0.01
nsteps=10000
niter=7
path=[]
for i in range(niter):
    random.seed(i)
    steps=[array([0,0])]
    for j in range(nsteps):
        th=pi/4*int(8*random.rand())
        steps.append(steps[-1]+stepsize*(-1 if random.rand() < 0.5 else 1)*array([cos(th),sin(th)]))
    path.append(steps)

path=array(path)
print 'average distance',sum([sqrt(p[-1,0]**2+p[-1,1]**2) for p in path])/niter

fig=figure()
ax=fig.add_subplot(111,aspect='equal')
for i in range(niter): plot(path[i,:,0],path[i,:,1])
xlim([-1,1])
ylim([-1,1])
savefig('2D_paths.pdf')
show()
"""
