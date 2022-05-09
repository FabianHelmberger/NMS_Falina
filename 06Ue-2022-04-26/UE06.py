#
from ast import Break
from cProfile import label
import numpy as np
import matplotlib.pyplot as plt

A = np.array([[6., -2., 3.], [1., 10., -4.], [-2., 1., 0.]])


def Gershgorin(A):
    dim = np.shape(A)[0]
    Radius = []

    out = []

    # loop over diagonal elemts
    for j in range(dim):
        tuple = []

        # loop oder cols and rows
        for dir in ['col', 'row']:
            if dir == 'row':
                A = A.T

            diff = 0
            for i in range(dim):

                # don't sum diagonals
                if i != j:
                    diff = diff + np.abs(A[j,i])

            tuple.append([diff])
        out.append(tuple)
    return out


sol = Gershgorin(A)
fig, axes = plt.subplots(figsize=(15, 15))

cls = ['red', 'green', 'blue']

for element, cl in zip(range(len(sol)), cls):
    for radius in sol[element]:
        
        diag = A[element, element]
        circle = plt.Circle((diag, 0), radius[0], alpha=0.3, color = cl)
        axes.add_patch(circle)
eig = np.linalg.eigvals(A)


y = np.zeros(3)
plt.scatter(eig, y)
plt.xlim([-20, 20])
plt.ylim([-20, 20])
plt.show()


# condition number 

# find the max / min center of center 
max_diag = np.max(np.diag(A))
min_diag = np.min(np.diag(A))

# find the corresponding indices
max_index = np.where(np.diag(A) == max_diag)[0][0]
min_index = np.where(np.diag(A) == min_diag)[0][0]

# find the max radii at these centeres
max_radius = np.max(sol[max_index])
min_radius = np.min(sol[min_index])

# find the corresponding
max_lambda = max_diag + max_radius
min_lambda = min_diag - min_radius

cond = np.abs(max_lambda / min_lambda)
print(cond)





# 2nd Problem: power method

# define z0
z0 = np.array([1,1,1])

# Calculate the unit vector
y0 = z0 / np.linalg.norm(z0)

# Calculate the new vector
z1 = A@y0

# Calculate the first estimate of Lambda
lambda_ = np.inner(y0, z1) / (np.inner(y0, y0))

# Define a convergence criteria
tol = 1e-6
err = 1 / tol
max_iter = 100
L = [lambda_]
E = []
k = 0

z = z1
while (err > tol and k < max_iter):
    y = z / np.linalg.norm(z)
    z = A@y
    lambda_old = lambda_
    lambda_ = np.inner(y, z) / np.inner(y, y)

    err = np.abs(lambda_old - lambda_)
    print(err)
    k = k+1
    L.append(lambda_)
    E.append(err)
print("Convergence (", err, ") was reached after ", k, " Iterations")
plt.scatter(range(np.shape(L)[0]), L)
plt.scatter(range(np.shape(E)[0]), E)
plt.show()



# 2nd Problem: QR Algo
def qr_step(A, q):
    return q.T@A@q

tol = 1e-6
err = 1/tol
max_iter = 100
k = 0
eigenvalues = []

while (err > tol and k < max_iter):
    q = np.linalg.qr(A)[0]
    A_ = A
    A = qr_step(A, q)
    err = np.abs(np.max(A-A_))
    eigenvalues.append(np.diagonal(A))
    k = k+1
    
print(k, ' Steps are needed to reach error = ', tol)

eigenvalues = np.array(eigenvalues)
plt.plot(eigenvalues[:, 0], label=r'$\lambda_1$')
plt.plot(eigenvalues[:, 1], label=r'$\lambda_2$')
plt.plot(eigenvalues[:, 2], label=r'$\lambda_3$')
plt.xlabel('iteration')
plt.ylabel('Eigenvalue')
plt.legend()
plt.show()