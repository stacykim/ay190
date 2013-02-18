# advect.py
# created by stacy kim
#
# Contains a suite of solvers for the advection equation. The solution is not 
# computed at the boundary points.


def upwind(y):
    """Implements the upwind method (O[h] in space and time)."""
    return [y[i]-v*dt/dx*(y[i]-y[i-1]) for i in range(1,len(y)-1)]

def downwind(y):
    """Implements the downwind method (O[h] in space and time)."""
    return [y[i]-v*dt/dx*(y[i+1]-y[i]) for i in range(1,len(y)-1)]

def ftcs(y):
    """Implements the FTCS method (O[h^2] in space, O[h] in time)."""
    return [y[i]-v*dt/(2*dx)*(y[i+1]-y[i-1]) for i in range(1,len(y)-1)]

def lax_friedrich(y):
    """Implements the Lax-Friedrich method (O[h^2] in space, O[h] in time)."""
    return [0.5*(y[i+1]+y[i-1]) - v*dt/(2*dx)*(y[i+1]-y[i-1]) 
            for i in range(1,len(y)-1)]

def leapfrog(y,yold2):
    """Implements the leapfrog method (O[h^2] in space and time)."""
    return [yold2[i] - v*dt/dx*(y[i+1]-y[i-1]) for i in range(1,len(y)-1)]

def lax_wendroff(y):
    """Implements the Lax-Wendroff method (O[h^2] in space and time)."""
    return [y[i]-v*dt/(2*dx)*(y[i+1]-y[i-1])+0.5*(v*dt/dx)**2*(y[i-1]-2*y[i]+y[i+1])
            for i in range(1,len(y)-1)]
