# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:32:00 2024

@author: aelxa
"""
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from numpy.linalg import norm 
from newton_and_quasinewton_script import newton_method_nd
from newton_and_quasinewton_script import lazy_newton_method_nd
from newton_and_quasinewton_script import broyden_method_nd
from interp import dividedDiffTable
from interp import evalDDpoly
from interp import eval_lagrange


matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.dpi'] = 300

def evalF(x): 

    F = np.zeros(2)
    
    F[0] = x[0]**2 + x[1]**2 - 4
    F[1] = np.exp(x[0]) + x[1] - 1
    
    return F 

def evalJ(x): 

    J = np.array([[2*x[0], 2*x[1]], 
       [np.exp(x[0]), 1]])
    return J

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


#1a
x0 = np.array([1.,1.])
tol = 1e-8J
Nmax = 100
[xstar1a, xlist1a, ier1a, its1a] = LazyNewton(x0, tol, Nmax)
print('Lazy Newton converged to ', xstar1a, ' in ', its1a, 'iterations')

#1c
[xstar1c, xlist1c,ier1c, its1c] = Newton(x0, tol, Nmax)
print('Newton converged to ', xstar1c, ' in ', its1c, 'iterations')

def error(its, xlist, xstar):
    err = np.zeros((its+1, 2))
    for i in range(its+1):
        err[i] = norm(xlist[i] - xstar)
    return err
    
err_ln = error(its1a, xlist1a, xstar1a)
err_n = error(its1c, xlist1c, xstar1c)

plt.figure()
plt.title('Lazy Newton error vs iteration')
pltits_ln = np.arange(0,its1a+1)    
plt.xlabel('Iteration')
plt.ylabel('log_10(norm(error))')
plt.semilogy(pltits_ln, err_ln)
plt.show()

plt.figure()
plt.title('Newton error vs iteration')
pltits_n = np.arange(0,its1c+1)    
plt.xlabel('Iteration')
plt.ylabel('log_10(norm(error))')
plt.semilogy(pltits_n, err_n)
plt.show()



#2
def F(x):
    return np.array([x[0]**2 + x[1]**2 - 4 ,  np.exp(x[0]) + x[1] - 1]);

def JF(x):
    return np.array([[2*x[0], 2*x[1]], 
       [np.exp(x[0]), 1]]);

# Apply Newton Method:
x0i = np.array([1.0,1.0]); tol=1e-14; nmax=100;
(rN,rnN,nfN,nJN) = newton_method_nd(F,JF,x0i,tol,nmax,True);
print(rN)

# Apply Lazy Newton (chord iteration)
nmax=1000;
(rLN,rnLN,nfLN,nJLN) = lazy_newton_method_nd(F,JF,x0i,tol,nmax,True);

# Apply Broyden Method
Bmat='fwd'; B0 = JF(x0); nmax=100;
(rB,rnB,nfB) = broyden_method_nd(F,B0,x0i,tol,nmax,Bmat,True);

# Plots and comparisons
numN = rnN.shape[0];
errN = np.max(np.abs(rnN[0:(numN-1)]-rN),1);
plt.plot(np.arange(numN-1),np.log10(errN+1e-18),'b-',label='Newton');
plt.title('Newton iteration log10|r-rn|');
plt.legend();
plt.show();

numB = rnB.shape[0];
errB = np.max(np.abs(rnB[0:(numB-1)]-rN),1);
plt.plot(np.arange(numN-1),np.log10(errN+1e-18),'b-',label='Newton');
plt.plot(np.arange(numB-1),np.log10(errB+1e-18),'g-',label='Broyden');
plt.title('Newton and Broyden iterations log10|r-rn|');
plt.legend();
plt.show();
 
#Lazy Newton does not work here, so the plotting of it is excluded here

#ii
x0ii = np.array([1.0, -1.0])

# Apply Newton Method:
x0 = x0ii; tol=1e-14; nmax=100;
(rN,rnN,nfN,nJN) = newton_method_nd(F,JF,x0ii,tol,nmax,True);
print(rN)

# Apply Lazy Newton (chord iteration)
nmax=1000;
(rLN,rnLN,nfLN,nJLN) = lazy_newton_method_nd(F,JF,x0ii,tol,nmax,True);

# Apply Broyden Method
Bmat='fwd'; B0 = JF(x0); nmax=100;
(rB,rnB,nfB) = broyden_method_nd(F,B0,x0ii,tol,nmax,Bmat,True);

#Plots and comparisons
numN = rnN.shape[0];
errN = np.max(np.abs(rnN[0:(numN-1)]-rN),1);
plt.plot(np.arange(numN-1),np.log10(errN+1e-18),'b-',label='Newton');
plt.title('Newton iteration log10|r-rn|');
plt.legend();
plt.show();

numB = rnB.shape[0];
errB = np.max(np.abs(rnB[0:(numB-1)]-rN),1);
plt.plot(np.arange(numN-1),np.log10(errN+1e-18),'b-',label='Newton');
plt.plot(np.arange(numB-1),np.log10(errB+1e-18),'g-',label='Broyden');
plt.title('Newton and Broyden iterations log10|r-rn|');
plt.legend();
plt.show();

numLN = rnLN.shape[0];
errLN = np.max(np.abs(rnLN[0:(numLN-1)]-rN),1);
plt.plot(np.arange(numN-1),np.log10(errN+1e-18),'b-',label='Newton');
plt.plot(np.arange(numB-1),np.log10(errB+1e-18),'g-',label='Broyden');
plt.plot(np.arange(numLN-1),np.log10(errLN+1e-18),'r-',label='Lazy Newton');
plt.title('Newton, Broyden and Lazy Newton iterations log10|r-rn|');
plt.legend();
plt.show();

#iii
x0iii = np.array([0.0, 0.0])

# Apply Newton Method:
x0 = x0iii; tol=1e-14; nmax=100;
(rN,rnN,nfN,nJN) = newton_method_nd(F,JF,x0iii,tol,nmax,True);
print(rN)

# Apply Lazy Newton (chord iteration)
nmax=1000;
(rLN,rnLN,nfLN,nJLN) = lazy_newton_method_nd(F,JF,x0iii,tol,nmax,True);

# Apply Broyden Method
Bmat='fwd'; B0 = JF(x0); nmax=100;
(rB,rnB,nfB) = broyden_method_nd(F,B0,x0iii,tol,nmax,Bmat,True);

# Plots and comparisons
numN = rnN.shape[0];
errN = np.max(np.abs(rnN[0:(numN-1)]-rN),1);
plt.plot(np.arange(numN-1),np.log10(errN+1e-18),'b-o',label='Newton');
plt.title('Newton iteration log10|r-rn|');
plt.legend();
plt.show();

numB = rnB.shape[0];
errB = np.max(np.abs(rnB[0:(numB-1)]-rN),1);
plt.plot(np.arange(numN-1),np.log10(errN+1e-18),'b-',label='Newton');
plt.plot(np.arange(numB-1),np.log10(errB+1e-18),'g-',label='Broyden');
plt.title('Newton and Broyden iterations log10|r-rn|');
plt.legend();
plt.show();

numLN = rnLN.shape[0];
errLN = np.max(np.abs(rnLN[0:(numLN-1)]-rN),1);
plt.plot(np.arange(numN-1),np.log10(errN+1e-18),'b-',label='Newton');
plt.plot(np.arange(numB-1),np.log10(errB+1e-18),'g-',label='Broyden');
plt.plot(np.arange(numLN-1),np.log10(errLN+1e-18),'r-',label='Lazy Newton');
plt.title('Newton, Broyden and Lazy Newton iterations log10|r-rn|');
plt.legend();
plt.show();

#4
N = 4
''' interval'''
a = 0
b = 8


''' create equispaced interpolation nodes'''
xint = [0,2,3,5,8]

''' create interpolation data'''
yint = [-125, -27, -8, 0, 27]

''' create points for evaluating the Lagrange interpolating polynomial'''
Neval = 1000
xeval = np.linspace(a,b,Neval+1)
yeval_l= np.zeros(Neval+1)
yeval_dd = np.zeros(Neval+1)

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
      
fex = (xeval-5)**3   

plt.figure()    
plt.plot(xeval,fex,'r-')
plt.plot(xeval,yeval_l,'b--') 
plt.plot(xeval,yeval_dd,'c--')
plt.legend()

plt.figure() 
err_l = abs(yeval_l-fex)
err_dd = abs(yeval_dd-fex)
plt.semilogy(xeval,err_l,'r--',label='lagrange')
plt.semilogy(xeval,err_dd,'b--',label='Newton DD')
plt.legend()
plt.show()