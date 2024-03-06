# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 10:32:30 2024

@author: aelxa
"""
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from interp import evalDDpoly
from interp import eval_lagrange
from interp import dividedDiffTable
import matplotlib

matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.dpi'] = 300

def monomial(xint,yint):
    V = np.vander(xint)
    a = np.linalg.solve(V,yint)
    return a


def driver():

    f = lambda x: np.sinc(5*x)

    N = 6

    
    #N=18 is where p(x) >100
    ''' interval'''
    a = -1
    b = 1
   
   
    ''' create equispaced interpolation nodes'''
    xint = np.linspace(a,b,N+1)
    
    ''' create interpolation data'''
    yint = f(xint)
    
    Neval = 1000
    xeval = np.linspace(a,b,Neval+1)
    yeval_l= np.zeros(Neval+1)
    yeval_dd = np.zeros(Neval+1)
    yeval_monom = np.zeros(Neval+1)
    
    
    a = monomial(xint, yint)
    yeval_monom = np.polyval(a, xeval)
    
    
   
  
    '''Initialize and populate the first columns of the 
     divided difference matrix. We will pass the x vector'''
    y = np.zeros( (N+1, N+1) )
     
    for j in range(N+1):
       y[j][0]  = yint[j]

    y = dividedDiffTable(xint, y, N+1)
    ''' evaluate lagrange poly '''
    for kk in range(Neval+1):
       yeval_l[kk] = eval_lagrange(xeval[kk],xint,yint,N)
       yeval_dd[kk] = evalDDpoly(xeval[kk],xint,y,N)
          

    ''' create vector with exact values'''
    fex = f(xeval)
 
    print(max(yeval_l))
    print(max(yeval_monom))
    print(max(yeval_dd))

    plt.figure() 
    plt.title('N=6: f(x) vs p(x)')
    plt.plot(xeval,fex,'r-', label='actual')
    plt.plot(xeval,yeval_l,'b-.', label='lagrange') 
    plt.plot(xeval,yeval_monom, 'k--', label = 'Monomial')
    plt.plot(xeval,yeval_dd,'g:', label='Newton DD')
    plt.legend(loc = "lower right")
    plt.show()

    plt.figure() 
    plt.title('N=6: error of p(x)')
    err_l = abs(yeval_l-fex)
    err_dd = abs(yeval_dd-fex)
    err_m = abs(yeval_monom - fex)
    plt.semilogy(xeval,err_l,'b-.',label='lagrange')
    plt.semilogy(xeval,err_m,'k--',label='Monomial')
    plt.semilogy(xeval,err_dd,'g:',label='Newton DD')
    plt.legend(loc = "lower left")
    plt.show()



driver()


