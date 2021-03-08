# This little code fragment implements a fourth order compact finite difference scheme with edge formulae
# RMC May 2019 U. Alcala, in Alcala de Henares

from numpy import linalg as LA
import numpy as np # Already done earlier in the notebook
from scipy import sparse
import scipy as sp # Already done earlier in the notebook

n = 673
if n<= 3:
    raise Exception('Sorry, friend, but this code needs 4 or more pieces of data')

x = np.array([np.cos(np.pi*(n-1-i)/(n-1)) for i in range(n)])
h = np.diff(x)

y = np.zeros(x.shape)

# Build tridiagonal matrix (floating point type, by contagion)
# subs == vector of subdiagonal entries
# diags == vector on diagonal entries
# supers == vector of superdiagonal entries
# Matrix as in the paper Compact Finite Differences and Cubic Splines
diags = np.array([4.0 for i in range(n)])
diags[0] = h[1]*(h[1]+h[2])/(h[0]+h[1])/(h[0]+h[1]+h[2])
diags[-1] = h[-2]*(h[-3]+h[-2])/(h[-2]+h[-1])/(h[-3]+h[-2]+h[-1])

subs = np.array([4.0*h[i+1]**2/(h[i]+h[i+1])**2 for i in range(n-2)]+[1])

supers = np.array([1]+[4.0*h[i]**2/(h[i]+h[i+1])**2 for i in range(n-2)])

# This factors the tridiagonal matrix.  We need only adjust the
# lower vector and the diagonal vector (the matrix is positive definite)
ell = np.zeros(subs.shape)
de  = np.zeros(diags.shape)
de[0] = diags[0]
for i in range(1,n):
    ell[i-1] = subs[i-1]/de[i-1]
    de[i]    = diags[i] - ell[i-1]*supers[i-1]

# Some made-up function to differentiate (to test)
fndata = np.zeros(diags.shape)
def fn(x):
    return 1.0/(1+x*x) # np.cos(x)
def dfn(x):
    return -2*x/(1+x*x)**2 #-1*np.sin(x)

for i in range(n):
    fndata[i] = fn(x[i])

# We will carry out the B.b operation manually; this will save storing the vectors
b = np.zeros(fndata.shape)

# Hand-translated from the Matlab (quite surprising that it worked, really)
a0 = (4*h[0]**2+6*h[0]*h[1]+3*h[0]*h[2]+2*h[1]**2+2*h[1]*h[2])*h[1]*(h[2]+h[1])/(h[0]+h[1])**2/(h[0]+h[1]+h[2])**2/h[0]
a1 = 1/h[0]*((-2*h[1]+h[0])*h[2]+2*h[1]*(-h[1]+h[0]))/h[1]/(h[1]+h[2])
a2 = -h[0]**2*(h[2]+h[1])/(h[1]+h[0])**2/(h[1])/h[2]
a3 = h[0]**2*h[1]/(h[2]+h[1]+h[0])**2/(h[2]+h[1])/h[2]
b[0] = -(a0*fndata[0] + a1*fndata[1] + a2*fndata[2] + a3*fndata[3])
for i in range(1,n-1):
    a0 = 4*h[i]**2*2*(2*h[i-1]+h[i])/h[i-1]/(h[i-1]+h[i])**3
    a1 = -8*(-h[i-1]+h[i])/h[i-1]/h[i]
    a2 = -8*h[i-1]**2*(h[i-1]+2*h[i])/(h[i-1]+h[i])**3/h[i]
    b[i] = -(a0*fndata[i-1] + a1*fndata[i] + a2*fndata[i+1])
    
a0 = -h[-2]*h[-1]**2/(h[-3]+h[-2]+h[-1])**2/(h[-3]+h[-2])/h[-3]
a1 = h[-1]**2/(h[-2]+h[-1])**2/h[-2]/h[-3]*(h[-2]+h[-3])
a2 = 1/h[-1]*((2*h[-2]-h[-1])*h[-3]+2*h[-2]*(h[-2]-h[-1]))/h[-2]/(h[-3]+h[-2])
a3 = -(4*h[-1]**2+(6*h[-2]+3*h[-3])*h[-1]+2*h[-2]*(h[-3]+h[-2]))/h[-1]/(h[-2]+h[-1])**2/(h[-3]+h[-2]+h[-1])**2*(h[-2])*(h[-2]+h[-3])

b[-1] = -(a0*fndata[-4] + a1*fndata[-3] + a2*fndata[-2] + a3*fndata[-1])

# This solves the tridiagonal system, given the factoring above.
# I should really peel those out into TRIDECOMP and TRISOLVE,
# implementing the so-called Thomas algorithm in a reasonably
# maintainable and re-usable way. Well, this is at this point
# just an exercise, and who will ever see this code except me?
#
# Oh, wait.  That never happens.  This code is going to escape into the wild,
# isn't it.
#
# Oops. nvm
#
y = np.zeros(b.shape)
y[0] = b[0]
for i in range(1,n):
    y[i] = b[i] - ell[i-1]*y[i-1]

sol = np.zeros(x.shape)
sol[-1] = y[-1]/de[-1]

for i in range(n-2,-1,-1): # start, stop, step (goofy language)
    sol[i] = (y[i]-supers[i]*sol[i+1])/de[i]

dy = np.zeros(x.shape)
for i in range(n):
    dy[i] = dfn(x[i])

print(LA.norm(sol-dy,2)*(n/2)**4) # This is the scaled error of derivatives at all the nodes