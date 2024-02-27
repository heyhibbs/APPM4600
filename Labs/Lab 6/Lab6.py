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
fforward = (np.cos((np.pi/2)+h)-np.cos(np.pi/2))/h
fcentered = (np.cos((np.pi/2)+h)-np.cos((np.pi/2)-h))/(2*h)

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
    
    F[0] = 4*x[0]**2+x[1]**2-4
    F[1] = x[0]+x[1] - np.sin(x[0]-x[1])
    return F
    
def evalJ(x): 

    J = np.array([[8*x[0], 2*x[1]], 
                  [1-np.cos(x[0]-x[1]), 1+np.cos(x[0]-x[1])]])
    return J


def SlackerNewton(x0,tol,Nmax):

    ''' Slacker Newton = recompute the Jacobian given a condition: Here we will use every 3 iterations '''
    ''' inputs: x0 = initial guess, tol = tolerance, Nmax = max its'''
    ''' Outputs: xstar= approx root, ier = error message, its = num its'''

    J = evalJ(x0)
    Jinv = inv(J)
    for its in range(Nmax):
        
       F = evalF(x0)
       x1 = x0 - Jinv.dot(F)
       
       if (its-1)%3 == 0:
           J = evalJ(x1)
           Jinv = inv(J)
           F = evalF(x1)
      
       if (norm(x1-x0) < tol):
           xstar = x1
           ier =0
           return[xstar, ier,its]
           
       x0 = x1
    
    xstar = x1
    ier = 1
    return[xstar,ier,its]   

# use routines
x0 = np.array([1,0])

Nmax = 100
tol = 1e-10

t = time.time()
for j in range(20):
  [xstar,ier,its] =  SlackerNewton(x0,tol,Nmax)
elapsed = time.time()-t
print(xstar)