# integrate.py
# created 10/5/11 by stacy kim
#
# A suite of quadrature routines, including
#   Newton-Cotes methods (piecewise polynomial interpolation):
#     midpoint rule
#     trapezoidal rule
#     Simpson's rule
#   Gaussian quadrature
#     Gauss-Legendre W(x) = 1
#     Gauss-Laguerre W(x) = x^c * exp(-x)  # here c=0
# All routines integrate the function f over the interval [a,b] over N equally
# spaced sub-intervals.
#
# modified 01/18-20/13
#   made more concise
#   implemented gaussian quadrature routines
#   added new routine integrate_to_accuracy()


import math
import numpy as np
import scipy.special as sp



# NEWTON-COTES METHODS ------------------------------------------------------------

def midpoint(f, a, b, N):
    """Integrates the given function via the midpoint rule."""
    h=float(b-a)/N
    return h*sum([f(x) for x in np.arange(a+h/2,b-h/2,h)])


def trapezoid(f, a, b, N):
    """Integrates the given function via the trapezoidal rule."""
    h=float(b-a)/N
    return h*(sum([f(x) for x in np.arange(a,b,h)-0.5*(f(a)+f(b))]))


def simpson(f, a, b, N):
    """Approximates the integral of a given function via Simpson's formula."""
    h=float(b-a)/N
    return h/6*sum([f(x)+4*f(x+h/2)+f(x) for x in np.arange(a,b,h)])



# GAUSSIAN QUADRATURE METHODS -----------------------------------------------------

def gauss_legendre(f, a, b, n):
    [roots,weights] = sp.p_roots(n,0)
    return (b-a)/2.*sum([weights[i]*f((b-a)/2.*roots[i]+(a+b)/2.) for i in range(n)])


def gauss_laguerre(f, a, b, n):
    if a != 0 and b != float('int'):
        raise ValueError('Gauss-Laguerre Quadrature integrates over [a,0).')

    w = lambda x: math.exp(x)  # inverse weighting function
    [roots,weights] = sp.l_roots(n,0)
    return sum([weights[i]*w(roots[i])*f(roots[i]) for i in range(n)])



# INTEGRATE TO ACCURACY -----------------------------------------------------------

N=2 # static variable
def integrate_to_accuracy(method, f, a, b, accuracy):
    """Integrates the given function by the given quadrature routine to the given 
    accuracy.  Returns the approximate value of the integral and the number of
    sub-intervals or nodes used to obtain the given accuracy."""
    I2=method(f,a,b,N)
    error=9e99

    while (error>accuracy):
        N*=2
        I1=I2
        I2=method(f,a,b,N)
        error=abs((I1-I2)/I1)

    #print "converged to accuracy",accuracy,"with",N,"subintervals"
    return [I2,N]
