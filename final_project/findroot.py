# findroot.py
# created 11/17/11 by stacy kim
#
# Contains a suite of routines to find roots that implement the methods
#    Newton-Raphson
#    secant
#    bisection.
#
#
# modified 1/26/13: removed support to output to file, instead returning array
# modified 1/28/13 to return number of iterations required to find solution
# modified 1/29/13: added method to find all rational roots to a polynomial


import sys
import math


# FILE OUTPUT ROUTINES -----------------------------------------------------------

files=None

def set_output(h,c,fns):
    """
    Enables output to files and sets file names for all 3 methods to those
    given in the array fns, in the order bisection, newton_raphson, secant.
    """
    global files

    # Open output files and write headers
    files=[0]*3

    for i in range(3):
        files[i]=open(fns[i],'w')
        files[i].write(fns[i])
        files[i].write('\nh={0},c={1},err=1e-4\n{2}{3}\n'.format(h,c,'x'.rjust(19),'precision'.rjust(19)))


def close_files():
    """Closes output files, if any."""
    if (files==None): return

    print "Finished writing output files."
    for i in range(3):
        files[i].close()


# ROOT-FINDING ROUTINES --------------------------------------------------------

def bisection(x1,x2,f,err):

    if (f(x1)*f(x2) > 0):
        raise ValueError('f(x1) and f(x2) must have opposite signs!\n'
                         '   f({0})={1},  f({2})={3}'.format(x1,f(x1),x2,f(x2)))

    if (x2==min(x1,x2)): x1,x2 = x2,x1

    count=0
    x=(x1+x2)/2.0
    while(abs(f(x)) > err):
        if (f(x)*f(x1)>0): x1=x
        else: x2=x
        x=(x1+x2)/2.0

        count+=1
        if (count>=990):
            print 'f({0})={1}, f({2})={3}'.format(x1,f(x1),x2,f(x2))
            if (count==1000):
                print 'Failed to converge after 1000 iterations.'
                sys.exit()

    return x,count


def newton_raphson(x,f,df,err):
    count=0
    while(abs(f(x)) > err): 
        x=x-f(x)/df(x)
        count+=1
    return x,abs(f(x)/df(x)),count


def secant(x2,x,f,err):
    count=0
    while(abs(f(x)) > err):
        x1=x2
        x2=x
        x=x2-f(x2)*(x2-x1)/(f(x2)-f(x1))
        count+=1
        if count % 10000 == 0: print count,':',x1,x2,x,f(x)
        
    df=(f(x2)-f(x))/(x2-x)
    return x,abs(f(x)/df),count


def find_all_roots(coeff0):
    """
    Computes all integer or rational roots of the given function 
    f = sum([coeff[i]*x**i for i=0..len(coeff)]).
    """
    coeff=coeff0

    roots=[]
    while len(roots) < len(coeff0)-1:
        # find zero of polynomial
        f  = lambda x: sum([coeff[i] * x**i for i in range(len(coeff))])
        root,err,niter=secant(0.0,1.0,f,1e-15)
        roots.append(root)

        #synthetic division
        a=coeff[-1]
        for i in range(len(coeff)-2,-1,-1): coeff[i]=a=coeff[i]+a*root
        if not (abs(coeff[0]) <= 1e-10):
            print 'found false root x =',root,'! (',abs(coeff[0]),')'
            sys.exit()
        coeff=coeff[1:]

    return roots
