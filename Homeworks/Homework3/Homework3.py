# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:16:49 2024

@author: aelxa

Homework 3-APPM4600
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axis import Axis
from bisection_example import bisection
import matplotlib
from fixedpt_example import fixedpt

matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.dpi'] = 300

#1c
f = lambda x: 2*x-1-np.sin(x)
[astar,ier, count]  = bisection(f, -np.pi, np.pi, 1e-8)
print('the approximate root is',astar)
print('the error message reads:',ier)
print('f(astar) =', f(astar))
print('this method took', count, 'iterations')

#2a
f2 = lambda x: (x-5)**9
a = 4.82
b = 5.2
tol = 1e-4
[astar, ier, count] = bisection(f2,a,b,tol)
print('the approximate root is',astar)
print('the error message reads:',ier)
print('f(astar) =', f2(astar))
print('this method took', count, 'iterations')

#2b
f2expanded = lambda x: x**9 - 45*x**8 + 900*x**7 - 10500*x**6 + 78750*x**5 - 393750*x**4 + 1312500*x**3 - 2812500*x**2 + 3515625*x - 1953125
[astar, ier, count] = bisection(f2expanded,a,b,tol)
print('the approximate root is',astar)
print('the error message reads:',ier)
print('f(astar) =', f2expanded(astar))
print('this method took', count, 'iterations')

#3b
f = lambda x: x**3 +x - 4
a = 1
b = 4
tol = 10**(-3)
[astar, ier, count] = bisection (f,a,b,tol)
print('the approximate root is',astar)
print('the error message reads:',ier)
print('f(astar) =', f(astar))
print('this method took', count, 'iterations')

#driver()

#5a
x = np.linspace(-10, 10, num = 150)
f5 = x - 4*np.sin(2*x)-3


fig, ax = plt.subplots()
ax.plot(x,f5)
Axis.set(ax, xlabel = 'x', ylabel = 'f(x)')
ax.grid()
plt.show()

#5b
f5b = lambda x: -np.sin(2*x) +5*x/4 - 3/4

x01 = 3
Nmax = 100
tol = 0.5*10**-10
[xstar1, ier] = fixedpt(f5b,x01,tol,Nmax)
print('Root 1 is', xstar1)

x02 = 0
[xstar2, ier] = fixedpt(f5b,x02,tol,Nmax)
print('Root 2 is', xstar2)

x03 = 1.5 #~1.732 cannot be found
[xstar3, ier] = fixedpt(f5b,x03,tol,Nmax)
print('Root 3 is', xstar3)

x04 = -1 #~-0.898 cannot be found
[xstar4, ier] = fixedpt(f5b,x04,tol,Nmax)
print('Root 4 is', xstar4)

x05 = 5 #~4.518 cannot be found
[xstar5, ier] = fixedpt(f5b,x05,tol,Nmax)
print('Root 5 is', xstar5)

x = np.linspace(-10, 10, num = 150)
df5dx = -2*np.cos(2*x) +5/4


fig, ax = plt.subplots()
plt.plot(x,f5, label = 'f(x) = x_n+1')
plt.plot(x,df5dx, label = 'f(x)' )
plt.xlabel('x')
plt.grid()
plt.legend()
plt.show()