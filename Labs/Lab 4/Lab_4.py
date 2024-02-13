# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 10:29:38 2024

@author: aelxa
"""
from fixedpt_example import fixedpt

def steffenson(f,a, tol, Nmax):
    count = 0
    p = [a]
    b = f(a)
    c = f(b)
    while count < Nmax:
        pn1 = a - (b-a)**2/(c-2*b+a)
        p.append(pn1)
        count = count + 1
        if (abs(p[count]-p[count-1]) <tol):
            return[count, p]

      
        
''' 
        while (count <Nmax):
           count = count +1
           x1 = f(x0)
           out.append(x1)
           if (abs(x1-x0) <tol):
              xstar = x1
              ier = 0
              return [xstar,ier,out]
           x0 = x1

        xstar = x1
        ier = 1
        return [xstar, ier, out]
'''