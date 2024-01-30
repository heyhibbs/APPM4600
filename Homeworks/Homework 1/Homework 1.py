# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 13:31:43 2024

@author: aelxa
"""
#Install necessary packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.dpi'] = 300

#1i.
intervals = int((2.08-1.92)/0.001)+1
x = np.linspace(1.920, 2.080, num = intervals)

p_coeff = x**9-18*x**8+144*x**7-672*x**6+2016*x**5 - 4032*x**4 + 5376*x**3 - 4608*x**2 + 2304 * x - 512

p_exp = (x-2)**9

plt.plot(x, p_coeff, label = 'via coefficients')
plt.plot(x, p_exp, label = 'via expression')
plt.legend(loc = 'lower right')
plt.xlabel('x')
plt.ylabel('p(x)')
plt.show()
#plt.savefig('Q1')

#5b
delta = np.logspace(-16,0,num = 17, base = 10)
x1 = np.pi
x2 = 10**6

y1 = np.cos(x1+delta)-np.cos(x1)
y2 = -2*np.sin(x1 + delta/2)*np.sin(delta/2)

y3 = np.cos(x2+delta)-np.cos(x2)
y4 = -2*np.sin(x2 + delta/2)*np.sin(delta/2)

plt.plot(delta, y1-y2, label = 'x=pi')
plt.plot(delta, y3-y4, label = 'x = 1000000')

plt.legend(loc='lower left')

plt.xscale('log')
plt.xlabel('delta')
plt.show()

#5c
taylor1 = (-delta*np.sin(x1))#-0.5*delta**2 
taylor2 = (-delta*np.sin(x2))#-0.5*delta**2 
plt.plot(delta, taylor1-y1, label = 'x = pi')
plt.plot(delta, taylor2-y2, label = 'x=1000000')
plt.xscale('log')
plt.xlabel('delta')
plt.ylabel('f(x) - Taylor Expansion')
plt.legend(loc = 'lower left')
plt.show()