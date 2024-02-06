# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 10:26:42 2024

@author: aelxa
"""
import numpy as np
from bisection_example import bisection

#1a
f = lambda x: x**2*(x-1)
bisection(f, 0.5, 2, 1e-5, 100)
#Root is ~1, no error


#1b
bisection(f, -1, 0.5, 1e-5, 100)
#Method cannot tell root

#1c
bisection(f, -1,2,1e-5,100)
#Root is ~1, no error

#For both a and c, the method converges to the root x = 1, however, because f(a)<0 and f(b)<0

#2a
f = lambda x: (x-1)*(x-3)*(x-5)
bisection(f, 0, 2.4, 1e-5, 100)
#Root is ~1

#2b
f = lambda x: (x-1)**1*(x-3)
bisection(f, 0, 2, 1e-5, 100)
#Root is ~1

#2c
f = lambda x: np.sin(x)
bisection(f, 0,0.1, 1e-5, 100)
