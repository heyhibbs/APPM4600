# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:06:50 2024

@author: aelxa
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
from scipy.special import gamma
import math

#2a
a = -5
b = 5
n = 50

def trapquad(a,b,f,n):
    t = np.linspace(a,b,n+1)
    h = (b-a)/n
    w = np.ones([n+1,1])
    w[0] = 0.5
    w[n] = 0.5
    w = w*h
    f = f(t)
    integrand = np.zeros([n+1,1])
    for i in range(0,n+1):
        integrand[i] = f[i]*w[i]
    return sum(integrand)
f = lambda s: 1/(1+s**2) 
trap_int = trapquad(a,b,f,n)
print('Trapezoidal rule evaluates to', trap_int)

def simpsonquad(a,b,f,n):
    t = np.linspace(a,b,n+1)
    h = (b-a)/n
    w = np.ones([n+1,1])
    for i in range(1,n):
        if i%2 == 0:
            w[i] = 2
        else:
            w[i] = 4
    w = w*h/3
    f = f(t)
    integrand = np.zeros([n+1,1])
    for i in range(0,n+1):
        integrand[i] = f[i]*w[i]
    return sum(integrand)
simp_int = simpsonquad(a,b,f,n)
print('Simpsons rule evaluates to', simp_int)

#2b
ddf = lambda s: 2*(-3*s**2+1)/(s**2+1)**3
def err_trap(a,b,ddf,n):
    t = np.linspace(a,b,1000)
    ddf = ddf(t)
    ddfalpha = max(abs(ddf))
    h = (b-a)
    return abs(ddfalpha*h**3/(12*n**2))
print('The maximum error for the trapezoidal rule is', err_trap(a,b,ddf,n))

d4f = lambda s: 24*(5*s**4-10*s**2+1)/(s**2+1)**5
def err_simp(a,b,d4f,n):
    t = np.linspace(a,b,1000)
    d4f = d4f(t)
    d4falpha = max(abs(d4f))
    h = (b-a)
    return abs(d4falpha*h**5/(180*n**4))
print('The maximum error for Simpsons rule is', err_simp(a,b,d4f,n))

actual = 2*np.arctan(5)
#Computed using Symbolab and symbolics

abs_err_trap = abs(trap_int-actual)
abs_err_simp = abs(simp_int - actual)

print('The error for trapezoidal rule is', abs_err_trap)
print('The error for Simpsons rule is', abs_err_simp)

#2c
sp_default, sp_def_err, sp_def_out = scipy.integrate.quad(f,a,b, epsabs = 10**(-6),full_output=True)
sp_104, sp_1014_err, sp_104_out = scipy.integrate.quad(f,a,b, epsabs = 10**(-4),full_output=True)
neval_def = sp_def_out['neval']
neval_104 = sp_104_out['neval']
print('Scipy quad with default tolerance evaluates to', sp_default)
print('Scipy quad with TOL < 10^-4 evaluates to', sp_104)
print('The number of evaluations for TOL = 10^-6 is', neval_def)
print('The number of evaluations for TOL = 10^-4 is', neval_104)

#4a
def trap_gamma(t):
    b = 5*t
    x = np.linspace(0,b)
    f = lambda x: x**(t-1.)*np.exp(-x)
    n = 175*t
    return trapquad(0,b,f,n)

for i in range(1,6):    
    print('Integral with trapezoidal method, t =', i*2, ':', trap_gamma(2*i))

for i in range(1,6): 
    print('Integral with scipy, t =', i*2, ':', gamma(2*i))


f = lambda x: x**(4-1.)
[x,w] = np.polynomial.laguerre.laggauss(3) #t=2
gammaGL = sum(f(x)*w)
print('Gauss-Laguerre evaluates to:', gammaGL)
print((gammaGL-gamma(4))/gamma(4))
    
