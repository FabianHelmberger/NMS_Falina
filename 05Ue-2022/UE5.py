from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
from time import process_time

# global parameters
n = 100
A = 10
B = 2*np.pi
phi_A = 1
phi_B = 0

# define charge distribution
def rho(x):
    return A*np.cos(B*x)

# exact sol
def phi(x):
    return A/B**2*np.cos(B*x) - x + 1 -A/B**2

# Thomas algorithm
def TA(off, diag, sol):
    n = len(sol)

    alpha = np.zeros(n)
    beta = np.zeros(n)
    out = np.zeros(n)
    out[0] = 1

    for i in range(n-1, 0, -1):
        alpha[i-1] = - off[i] / (diag[i] + off[i] * alpha[i])
        beta[i-1] =  (sol[i] - off[i] * beta[i])/(diag[i] + off[i] * alpha[i])

    for i in range(1, n):
        out[i] = alpha[i-1] * out[i-1] + beta[i-1]

    return out



def num_sol(n):
    # define x domain
    X = np.linspace(0, 1, num=n-1)
    h = X[1] - X[0]

    # define Matrix for 2nd deriv approx
    diag =  np.array([-2 for i in range(n)] * 1 / h ** 2)
    diag =  np.concatenate(([1], diag, [1]))

    off = np.array([1 for i in range(n-1)] * 1 / h ** 2)
    off =  np.concatenate(([0], off, [0]))

    F = -rho(X)
    F = np.concatenate(([1], F, [0]))
    
    num_sol = TA(off, diag, F)

    return num_sol


import time
# Plot the average acc as a func of n
avgacc  = []
avgtime = []
for n in range(10, 100):
    X = np.linspace(0, 1, num=n+1)

    t_start = time.time()
    num = num_sol(n)
    t_stop = time.time()

    exa = phi(X)
    avgacc.append(np.linalg.norm(num-exa) / n)
    avgtime.append(t_stop - t_start)

plt.plot(avgacc, label='average acc over n')
plt.plot(avgtime, label='average cpu-time over n')
plt.legend()
plt.xlabel('n')
plt.ylabel('acc, cpu-time')
plt.xscale("log")
plt.yscale("log")
plt.show()

##solutions
#t10_start = time.time()
#
#n=10
#exact_10=exact_sol=phi(np.linspace(0,1,n+1))
#num_10=num_sol(n)
#t10_stop = time.time()
#print("time elapsed for n=10:", t10_stop-t10_start) 
#
##function for difference of solutions
#diff_10=abs(exact_10-num_10)
#
#
#t100_start = time.time()
#n=100
#exact_100=exact_sol=phi(np.linspace(0,1,n+1))
#num_100=num_sol(n)
#t100_stop = time.time()
#print("time elapsed for n=100:", t100_stop-t100_start) 
#
##function for difference of solutions
#diff_100=abs(exact_100-num_100)
#