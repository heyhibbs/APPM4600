# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 20:19:33 2024

@author: aelxa
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.dpi'] = 300

#Pre-Lab
h = 0.01*2.0**(-np.arange(0,10))
fforward = (np.cos((np.pi/2)+h)-np.cos(np.pi/2))/h
fcentered = (np.cos((np.pi/2)+h)-np.cos((np.pi/2)-h))/(2*h)

fprime = -np.sin(np.pi/2)

eforward = abs((fforward-fprime)/fprime)
ecentered = abs((fcentered-fprime)/fprime)

plt.plot(h,fforward)
plt.xlabel('h')
plt.ylabel("f'")
plt.title("Forward Difference Approximation of f'")
plt.show()

plt.plot(h, fcentered)
plt.xlabel('h')
plt.ylabel("f'")
plt.title("Centered Difference Approximation of f'")
plt.show()

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

