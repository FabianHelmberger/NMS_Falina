import  numpy as np
import  matplotlib.pyplot as plt
import  time

from    cProfile    import label
from    scipy       import special


# global parameters
n = 100
A = 10
B = 2*np.pi

# define charge distribution
def rho(x):
    return A*np.cos(B*x)

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


# exact solution
def exa_sol(x):
    return A/B**2*np.cos(B*x) - x + 1 -A/B**2


# numerical solution
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



# Calculate the exact solution and plot the numerical solution against the
# exact solution for n = 10 and 100

n = 10
X = np.linspace(0, 1, num=n+1)
plt.xlabel('x')
plt.ylabel(r'$\phi$(x)')
plt.title('Potential $\phi$(x) over x')
plt.plot(X, num_sol(n),  label="numerical, n = 10")

n = 100
X = np.linspace(0, 1, num=n+1)
plt.plot(X, num_sol(n),  label="numerical, n = 100")
plt.plot(X, exa_sol(X),  label="exact")
plt.legend()
plt.show()


# Plot the average accuracy as a function of n
# Plot the average cpu-time as a function of n

# !! Here we use standard deviation as a measure for accuracy !!

# arrays holding std and cpu time
stds    = []
times   = []

for n in range(5, 5000):
    X = np.linspace(0, 1, num=n+1)

    t_start = time.time()
    num = num_sol(n)
    t_stop = time.time()
    exa = exa_sol(X)

    # compute standard deviation
    std_ = np.sqrt(1 / n * np.sum(np.square(num-exa)))

    # compute the cpu time
    time_ = t_stop - t_start

    stds.append(std_)
    times.append(time_)

plt.plot(stds,  label='sd' )
plt.plot(times, label='ct' )
plt.title('standard deviation (sd) & cpu time (ct) over n')
plt.legend()
plt.xlabel('n')
plt.ylabel('sd, ct')
plt.xscale("log")
plt.yscale("log")
plt.show()