import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rcParams['figure.dpi'] = 300

def driver():


    f = lambda x: 1/(1+(10*x)**2)

    N = 18;
    ''' interval'''
    a = -1;
    b = 1;

    ''' create equispaced interpolation nodes'''
    xint = np.linspace(a,b,N+1);

    ''' create interpolation data'''
    yint = f(xint);

    Neval = 1000;
    xeval = np.linspace(a,b,Neval+1);

    ''' Linear spline evaluation '''
    yeval_ls = eval_lin_spline(xeval,xint,yint,N);

    ''' create points for evaluating the Lagrange interpolating polynomial'''
    yeval_l= np.zeros(Neval+1)
    yeval_dd = np.zeros(Neval+1)

    '''Initialize and populate the first columns of the
     divided difference matrix. We will pass the x vector'''
    y = np.zeros( (N+1, N+1) )

    for j in range(N+1):
       y[j][0]  = yint[j]

    y = dividedDiffTable(xint, y, N+1)
    ''' evaluate lagrange poly '''
    for kk in range(Neval+1):
       yeval_l[kk] = eval_lagrange(xeval[kk],xint,yint,N)
    
       yeval_dd[kk] = evalDDpoly(xeval[kk],xint,y,N)

    ''' create vector with exact values'''
    fex = f(xeval)

    plt.figure()
    plt.plot(xeval,fex,'r-')
    plt.plot(xeval,yeval_l,'b--', label = 'Lagrange')
    plt.plot(xeval,yeval_dd,'c--', label = 'Newton DD')
    plt.plot(xeval,yeval_ls,'g--', label = 'Linear Spline')
    plt.legend()

    plt.figure()
    err_l = abs(yeval_l-fex)
    err_dd = abs(yeval_dd-fex)
    err_ls = abs(yeval_ls-fex)
    plt.semilogy(xeval,err_l,'r--',label='lagrange')
    plt.semilogy(xeval,err_dd,'b--',label='Newton DD')
    plt.semilogy(xeval,err_ls,'g--',label='lin spline')
    plt.legend()
    plt.show()

def eval_line(x,x0,y0,x1,y1):
    lin = (1/(x1-x0))*(y0*(x1-x) + y1*(x-x0));
    return lin;

def find_int(xeval,a,b):
    ind = np.where(np.logical_and(xeval>=a,xeval<=b));
    return ind;

def eval_lin_spline(xeval,xint,yint,N):
    Neval = len(xeval);
    yeval = np.zeros(Neval);

    for n in range(N):
        indn = find_int(xeval,xint[n],xint[n+1]);
        yeval[indn] = eval_line(xeval[indn],xint[n],yint[n],xint[n+1],yint[n+1]);

    return yeval;

def eval_lagrange(xeval,xint,yint,N):

    lj = np.ones(N+1)

    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.

    for jj in range(N+1):
       yeval = yeval + yint[jj]*lj[jj]

    return(yeval)


''' create divided difference matrix'''
def dividedDiffTable(x, y, n):

    for i in range(1, n):
        for j in range(n - i):
            y[j][i] = ((y[j][i - 1] - y[j + 1][i - 1]) /
                                     (x[j] - x[i + j]));
    return y;

def evalDDpoly(xval, xint,y,N):
    ''' evaluate the polynomial terms'''
    ptmp = np.zeros(N+1)

    ptmp[0] = 1.
    for j in range(N):
      ptmp[j+1] = ptmp[j]*(xval-xint[j])

    '''evaluate the divided difference polynomial'''
    yeval = 0.
    for j in range(N+1):
       yeval = yeval + y[0][j]*ptmp[j]

    return yeval



driver()
