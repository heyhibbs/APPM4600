# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 20:19:33 2024

@author: aelxa
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# get libraries
import numpy as np
import math
import time
from numpy.linalg import inv 
from numpy.linalg import norm 

matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.dpi'] = 300

#Pre-Lab
h = 0.01*2.0**(-np.arange(0,10))
f = lambda x: np.cos(x)
fforward = (f((np.pi/2)+h)-f(np.pi/2))/h

fcentered = (f((np.pi/2)+h)-f((np.pi/2)-h))/(2*h)

fprime = -np.sin(np.pi/2)

eforward = abs((fforward-fprime)/fprime)
ecentered = abs((fcentered-fprime)/fprime)

plt.loglog(1/h, eforward)
plt.ylabel('log(e)')
plt.xlabel('-log(h)')
plt.title('Forward Difference Error vs h')
plt.show()

plt.loglog(1/h, ecentered)
plt.ylabel('log(e)')
plt.xlabel('-log(h)')
plt.title('Centered Difference Error vs h')
plt.show()

orderf = (np.log(eforward[2])-np.log(eforward[1]))/(np.log(h[2])-np.log(h[1]))
orderc = (np.log(ecentered[2])-np.log(ecentered[1]))/(np.log(h[2])-np.log(h[1]))

#Lab 6
'''3.1: Recompute the Jacobian  '''
# define routines
def evalF(x): 

    F = np.zeros(2)
    
    F[0] = (x[0]-1)**2 + x[1]**2-1
    F[1] = (x[0]-2)**2 + (x[1]-1)**2-1
    return F
    
def evalJ(x): 

    J = np.array([[2*(x[0]-1), 2*x[1]], 
                  [2*(x[0]-2), 2*(x[1]-1)]])
    return J

def Newton(x0,tol,Nmax):

    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''
    xlist = np.zeros((Nmax+1,len(x0)));
    xlist[0] = x0;

    for its in range(Nmax):
       J = evalJ(x0);
       F = evalF(x0);

       x1 = x0 - np.linalg.solve(J,F);
       xlist[its+1]=x1;

       if (norm(x1-x0) < tol*norm(x0)):
           xstar = x1
           ier =0
           return[xstar, xlist,ier, its];

       x0 = x1

    xstar = x1
    ier = 1
    return[xstar,xlist,ier,its];


def LeaSlackerNewton(x0,tol,Nmax):

    ''' Slacker Newton = recompute the Jacobian given a condition: Here we will use every 3 iterations '''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    xlist = np.zeros((Nmax+1,len(x0)));
    xlist[0] = x0;

    J = evalJ(x0);
    for its in range(Nmax):
        if ((its+1)%3 == 0):
            J = evalJ(x0)
            
        F = evalF(x0)
        x1 = x0 - np.linalg.solve(J,F)
        xlist[its+1]=x1;

        if (norm(x1-x0) < tol*norm(x0)):
            xstar = x1
            ier =0
            return[xstar,xlist, ier,its]

        x0 = x1
    xstar = x1
    ier = 1
    return[xstar,xlist,ier,its];


def LazyNewton(x0,tol,Nmax):

    ''' Lazy Newton = use only the inverse of the Jacobian for initial guess'''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    xlist = np.zeros((Nmax+1,len(x0)));
    xlist[0] = x0;

    J = evalJ(x0);
    for its in range(Nmax):

       F = evalF(x0)
       x1 = x0 - np.linalg.solve(J,F);
       xlist[its+1]=x1;

       if (norm(x1-x0) < tol*norm(x0)):
           xstar = x1
           ier =0
           return[xstar,xlist, ier,its];

       x0 = x1

    xstar = x1
    ier = 1
    return[xstar,xlist,ier,its];


# use routines
x0 = np.array([1,2])

Nmax = 10
tol = 1e-10

[xstarNewton,xlistNewton, ierNewton,itsNewton] = Newton(x0, tol, Nmax)
[xstarSlack,xlistSlack, ierSlack,itsSlack] =  LeaSlackerNewton(x0,tol,Nmax)
[xstarLazy,xlistLazy, ierLazy,itsLazy] =  LazyNewton(x0,tol,Nmax)

print('Newton xstar', xstarNewton)
print('itsNewton', itsNewton)
print('Slacker xstar',xstarSlack)
print('itsSlack', itsSlack)
print('Lazy xstar',xstarLazy)
print('itsLazy', itsLazy)

