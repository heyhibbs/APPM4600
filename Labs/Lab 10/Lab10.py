# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 19:44:09 2024

@author: aelxa
"""
import numpy as np
import scipy as sc

#Prelab
def eval_legendre(n,x):
    phi = np.zeros(n+1)
    phi[0] = 1
    if n >= 1:
        phi[1] = x
    if n >= 2:
        for i in range(2,n+1):
            phi[i] = ((2*(i-1)+1)*x*phi[i-1]-(i-1)*phi[i-2])/((i-1)+1)
    if n == 0:
        p = phi[0]
    else:
        p = phi[0:n+1]
    return p
p = eval_legendre(4,2)
print(p)
