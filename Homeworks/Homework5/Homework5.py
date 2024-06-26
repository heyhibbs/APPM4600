# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 14:04:42 2024

@author: aelxa
HW 5
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math

matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.dpi'] = 300

#1

def eval_lagrange(xeval,xint,yint,N):

    lj = np.ones(N+1)
    #w = np.onoes(N+1)
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.
    
    for jj in range(N+1):
       yeval = yeval + yint[jj]*lj[jj]
  
    return(yeval)

def driver():


    f = lambda x: 1/(1+(16*x)**2)

    N = 4
    ''' interval'''
    a = -1
    b = 1
   
   
    ''' create equispaced interpolation nodes'''
    xint = np.linspace(a,b,N+1)
    
    ''' create interpolation data'''
    yint = f(xint)
    
    ''' create points for evaluating the Lagrange interpolating polynomial'''
    Neval = 1000
    xeval = np.linspace(a,b,Neval+1)
    yeval_l= np.zeros(Neval+1)
    #yeval_dd = np.zeros(Neval+1)
  
    '''Initialize and populate the first columns of the 
     divided difference matrix. We will pass the x vector'''
    y = np.zeros( (N+1, N+1) )
     
    for j in range(N+1):
       y[j][0]  = yint[j]
    ''' create vector with exact values'''
    fex = f(xeval)
       

    plt.figure()    
    plt.plot(xeval,fex,'ro-')
    plt.plot(xeval,yeval_l,'bs--') 
    #plt.plot(xeval,yeval_dd,'c.--')
    plt.legend()

    plt.figure() 
    err_l = abs(yeval_l-fex)
    #err_dd = abs(yeval_dd-fex)
    plt.semilogy(xeval,err_l,'ro--',label='lagrange')
    #plt.semilogy(xeval,err_dd,'bs--',label='Newton DD')
    plt.legend()
    plt.show()

driver()

#4
M = np.array([1, 1,1,1])
y = np.array([1,4,2,6])
W = np.zeros([4,4])
