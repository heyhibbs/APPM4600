import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
from time import perf_counter_ns


def driver():

     ''' create  matrix for testing different ways of solving a square 
     linear system'''

     '''' N = size of system'''
     N = 100
 
     ''' Right hand side'''
     
     b = np.random.rand(N,1)
     A = np.random.rand(N,N)
     
     start = perf_counter_ns()
     x = scila.solve(A,b)
     stop = perf_counter_ns()
     normaltime = stop-start
     print('Elapsed time for scipy solve:', normaltime, 'ns')
     
     startlu = perf_counter_ns()
     lu, piv = scila.lu_factor(A)
     stoplu = perf_counter_ns()
     lufactime = stoplu-startlu
     print('Elapsed time for LU factorization:', lufactime, 'ns')
     
     startlusolve = perf_counter_ns()
     xlu = scila.lu_solve((lu,piv),b)
     stoplusolve = perf_counter_ns()
     lusolvetime = stoplusolve-startlusolve
     print('Elapsed time for LU solve:', lusolvetime, 'ns')
     
     test = np.matmul(A,x)
     testlu = np.matmul(A,xlu)
     rlu = la.norm(testlu-b)
     r = la.norm(test-b)
     
     rhsnumb = (normaltime-lufactime)/lusolvetime
     if rhsnumb >=0:
         print('Number of rhs LU is faster:', rhsnumb)
     else:
         print('Number of rhs LU is faster:', 1)
     
     print(r)
     print(rlu)

     ''' Create an ill-conditioned rectangular matrix '''
     N = 10
     M = 5
     A = create_rect(N,M)     
     b = np.random.rand(N,1)


     
def create_rect(N,M):
     ''' this subroutine creates an ill-conditioned rectangular matrix'''
     a = np.linspace(1,10,M)
     d = 10**(-a)
     
     D2 = np.zeros((N,M))
     for j in range(0,M):
        D2[j,j] = d[j]
     
     '''' create matrices needed to manufacture the low rank matrix'''
     A = np.random.rand(N,N)
     Q1, R = la.qr(A)
     test = np.matmul(Q1,R)
     A =    np.random.rand(M,M)
     Q2,R = la.qr(A)
     test = np.matmul(Q2,R)
     
     B = np.matmul(Q1,D2)
     B = np.matmul(B,Q2)
     return B     
          
  
if __name__ == '__main__':
      # run the drivers only if this is called from the command line
      driver()       
