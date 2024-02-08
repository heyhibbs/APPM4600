# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 20:13:19 2024

@author: aelxa
"""

'''Homework 2 '''
import numpy as np
import math
import matplotlib.pyplot as  plt
import matplotlib
import random

matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.dpi'] = 300

#2b
A = np.matrix([[1, 1], [1+(10e-10), 1-(10e-10)]])
#print(A)
A = A*0.5
#print(A)
#A_inv = np.matrix([[1-10**10, 10**10], [1+(10**10), -(10**10)]])
A_inv = np.linalg.inv(A)

cond_num = np.linalg.cond(A, np.inf)
print('The condition number of A is', cond_num)

b = np.array([[1],[1]])
delta_b = np.array([[2e-5], [2e-5]])

b_tilda = b+delta_b

x = np.matmul(A_inv,b_tilda)
x_true = np.array([[1],[1]])

print(x)
rel_err = abs(x_true-x)/abs(x_true)
print("Relative error is", rel_err)

#3a

#3b
def f(x):
    y = math.exp(x)
    return y-1

y = f(9.999999995000000e-10)
print(y)

#3c
x = 9.999999995000000e-10
def f_taylor(x):
    y = 1 + x #+ (x**(2))/2 + (x**(3))/6 + (x**(4))/24 + x**5/120 + x**6/720 + x**7/5040 + x**8/math.factorial(8)
    return y-1
print(f_taylor(x))

rel_error = abs(10**(-9)-f_taylor(x))/abs(10**(-9))
print('rel_error is', rel_error)

#I'm losing my mind and I'll come back to this because it seems like the best I can get is like ~10^-8 digits of accuracy :(

#4a - there's no way this is not overcomplicated, but here we are
t = np.linspace(0,np.pi,30)
y = np.cos(t)

def suma(t,y):
    S = np.zeros(30)
    S[0] = t[0]*y[0]
    for i in range(1, 30):
        S[i] = t[i]*y[i]
    S = sum(S)
    return S
print(suma(t,y))

#4c
R = 1.2
deltar = 0.1
f = 15
p = 0
def x(theta, R, deltar, f, p):
    return R*(1+deltar*np.sin(f*theta+p))*np.cos(theta)

def y(theta,R, deltar, f, p):
    return R*(1+deltar*np.sin(f*theta+p))*np.sin(theta)

theta_range = np.linspace(0, 2*np.pi)
plt.plot(x(theta_range, R, deltar, f, p), y(theta_range,R, deltar, f, p))
plt.show()

#4b figure 2
for i in range(0,10):
    Ri = i
    deltari = 0.05
    fi = 2+i
    pi = random.uniform(0,2)
    plt.plot(x(theta_range, Ri, deltari, fi, pi), y(theta_range,Ri, deltari, fi, pi))
plt.show()
