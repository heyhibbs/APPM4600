# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 21:39:20 2024

@author: aelxa
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from numpy.linalg import inv

x = np.array ([0,2])
fx = ([2,4])

def construct_line(x,fx,alpha):
    m = (fx[1]-fx[0])/(x[1]-x[0])
    b =fx[0]-m*x[0]
    p = [m,b]
    leval = np.polyval(p,alpha)
    return leval

test = construct_line(x,fx,1)