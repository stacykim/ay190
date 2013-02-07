# rk4.py
# created 11/8/11 by stacy kim
#
# Computes a function given its derivative func with initial conditions x using the
# Runge-Kutta method.  2nd, 3rd, and 4th order routines and adaptive stepping are
# implemented.
#
#
# modified 11/10/11 to include adaptive step size control
# modified 1/24/12
#    corrected 2-step rungekutta call (t+h/2 in second call)
#    corrected error comparison in stepper routine (element vs. list comparison)
# modified 1/31/12
#    modified error comparison in stepper (max(list) vs. element)
#    removed calculation complete status output from driver
#    added debug print stmts to stepper to check error decrease
# modified 2/6/12: removed shebang line
# modified 2/12/12
#    added t0 to the driver's parameter list
#    create min step size (abort if reached)
# modified 1/21/13
#    renamed module from rk4 to rungekutta and 4th order routine
#    implemented 2nd and 3rd order runge-kutta methods
#    modified driver and stepper routines to use the rk method of the given order
#    switched support of vector arith. from homegrown Vector module to numpy


import sys, math

MAXITER=1000
MINH=1e-8
h=1.0


# TIME-STEPPING ROUTINES  -------------------------------------------------------

def driver(rk,x,t0,tf,h,f,err,fn):
    """
    Solves the differential equation 'f' using the given Runge-Kutta method 
    'rk' from t0 to tf.  Takes time steps of 'h' or adaptively steps to satisfy
    the given accuracy criteria 'err'.  Results are printed to the file fn.
    """

    ADAPTIVE = 1 if h==0 else 0

    # Open output file and write headers
    file=open(fn,'w')
    file.write(fn)
    file.write('\nh={0},t0={5},tf={1},err={4}\n{2}{3}'.format(h,tf,'t'.rjust(19),'h'.rjust(19),err,t0))
    for i in range(len(x)): file.write('x{0}'.format(i+1).rjust(19))
    file.write('\n')


    #Iteratively calculate position and velocities using given time steps
    t=t0
    i=0
    while (t<=tf):
        # Write results to file
        file.write('{0} {1} '.format(str(t).rjust(19),str(h).rjust(19)))
        for j in range(len(x)): file.write('{0} '.format(str(x[j]).rjust(19)))
        file.write('\n')

        # Calculate function value at next step
        if (ADAPTIVE):
            x,t,h=stepper(rk,x,t,f,err)
        else:
            t=h*i
            x=rk(x,t,h,f)
            i+=1

    file.close()

    return


def stepper(rk,x,t,f,err):
    """
    Computes the next step of the function f using the given Runge-Kutta method
    rk to the given accuracy 'err'.  Only called if adaptively time stepping.
    """

    global h
    h0=h
    n=0

    while (1):
        x1=rk(x,t,h,f)
        x2=rk(rk(x,t,h/2,f),t+h/2,h/2,f)
        
        if (max(abs(x2-x1)) > err[0]):
            h/=2
            if h < MINH:
                print 'Exceeded max precision.'
                sys.exit()
        else:
            h_old=h
            if max(abs(x2-x1)) < err[0]/10.: h*=2
            return x2,t+h_old,h_old

        n+=1
        if ((MAXITER-n)<10):
            print "abs(x2-x1)={0},err={1}\n".format(abs(x2-x1),err)
            if (n==MAXITER):
                print "Failed to converge within {0} iterations.".format(MAXITER)
                print "[x={0},err={1},h0={2},hf={3},t={4}]".format(x,abs(x2-x1),h0,h,t)
                sys.exit()



# RUNGE-KUTTA ROUTINES -----------------------------------------------------------

def rk4(x,t,h,f):
    """
    Computes the function value x with derivative f at the time t + h using
    the fourth order Runge-Kutta method.
    """

    k1=h*f(x,t)
    k2=h*f(x+k1/2,t+h/2)
    k3=h*f(x+k2/2,t+h/2)
    k4=h*f(x+k3,t+h)

    return x + (k1 + 2*k2 + 2*k3 + k4)/6


def rk3(x,t,h,f):
    """
    Computes the function value x with derivative f at the time t + h using
    the third order Runge-Kutta method.
    """
    
    k1=h*f(x,t)
    k2=h*f(x+k1/2,t+h/2)
    k3=h*f(x-k1+2*k2,t+h)

    return x + (k1 + 4*k2 + k3)/6


def rk2(x,t,h,f):
    """
    Computes the function value x with derivative f at the time t + h using
    the second order Runge-Kutta method.
    """

    k1=h*f(x,t)
    k2=h*f(x+k1/2,t+h/2)

    return x + k2
