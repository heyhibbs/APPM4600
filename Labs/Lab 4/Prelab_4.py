# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 21:53:04 2024

@author: aelxa
"""
import numpy as np
from fixedpt_example import fixedpt
import matplotlib.pyplot  as plt
import matplotlib

matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.dpi'] = 300

f1 = lambda x: (10/(x+4))**0.5
''' 
fixed point is alpha1 = 1.4987....
'''

Nmax = 100
tol = 1e-10
p =  1.3652300134140976

def order(out,p,alpha):
    lamb_duh = [abs(out[1]-out[0])/(abs(out[0]-p))**alpha]
    for i in range(1,(len(out)-1)):
        new = abs(out[i]-out[i-1])/(abs(out[i-1]-p))**alpha
        lamb_duh.append(new)
    return(lamb_duh)
    


x0 = 1.5
[xstar,ier,out] = fixedpt(f1,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f1(xstar):',f1(xstar))
print('Error message reads:',ier)

bama_0_5 = order(out,p,0.5)
bama1 = order(out, p ,1)
bama1_1 = order(out,p,1.5)
bama2 = order(out,p,2)

iterates = np.arange(0,len(out)-1)
plt.semilogy(iterates, bama_0_5, label='alpha = 0.5')
plt.semilogy(iterates, bama1, label='alpha = 1')
plt.semilogy(iterates, bama1_1, label='alpha = 1.5')
plt.semilogy(iterates, bama2, label='alpha = 2')
plt.xlabel('n')
plt.legend()
plt.show()

print('Number of iterates: ', len(out)+1)