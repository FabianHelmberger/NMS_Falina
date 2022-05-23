import numpy as np
import matplotlib.pyplot as plt
from sympy import E


################################################################################################
############################### NUMEROW RELATED FUNCTIONS ######################################
################################################################################################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Applying Numerow method to solve 
#!    d^2 Y(x)
#!    _______  =  - G(x)*Y(X) + S(X)
#!      d x^2
#!
#! by integration from inside ( X(1) ) to outside ( X(NR) )
#! ASSUMING PROPER INITIAL VALUES AT THE BOUNDARY!!!
#!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
def numerow_from_outside(y,g,s,n,dr):
#!###########################################################
#!   TODO: ENTER APPROPRIATE NUMEROW EXPRESSION HERE
#!###########################################################
    for i in np.arange(n-2, 1, -1):
        y[i+1]=(2*y[i]*(1-5*dr**2*g[i]/12)-y[i-1]*(1+dr**2*g[i-1]/12)+dr**2*(s[i+1]+10*s[i]+s[i-1])/12)/(1+dr**2*g[i+1]/12)
    return y

def numerow_from_inside(y,g,s,n,dr):
#!###########################################################
#!   TODO: ENTER APPROPRIATE NUMEROW EXPRESSION HERE
#!###########################################################
    for i in np.arange(2, n-1):
        y[i-1]=-(y[i+1]*(1+dr**2*g[i+1]/12)-2*y[i]*(1-5*dr**2*g[i]/12)-dr**2*(s[i+1]+10*s[i]+s[i-1])/12)/(1+dr**2*g[i-1]/12)
    return y


#######################################################################################
############################ Poisson-solver related functions  ########################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!
#! solve the  Poisson equation
#!
#!    d^2 UPOT(x)
#!    ___________  =  - URHO(X)
#!      d x^2
#!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#######################################################################################

def init_upot(r,n):
    upot=np.zeros(n)
#!###########################################################
#!   TODO: SET PROPEER INITIAL VALUES FOR THE POTENTIAL
#!         NEEDED BY NUMEROW METHOD (2 lines)
#!###########################################################
    upot[-1] = 1
    upot[-2] = 1
    return upot


def calc_urho(u,r,n):
    urho=np.zeros(n)
#!###########################################################
#!   TODO: CALCULATE CHARGE DENSITY HERE (2 lines)
#!###########################################################
    for i in range(n):
        urho[i] = u[i]**2 / r[i]
    return urho


def solve_poisson(upsi,erep,r,n,dr):
# arrays used for numerow method

    # array needed for numerow method
    g=np.zeros(n)

    upot=np.zeros(n)
#!###########################################################
#!   TODO: CALL THE THREE REQUIRED FUNCTIONS IN THE CORRECT ORDER
#!         TO OBTAIN THE SOLUTION FOR THE POISSON EQUATION HERE (3 lines)
#!###########################################################
    urho = calc_urho(upsi, r, n)
    upot = init_upot(r, n)
    y = numerow_from_outside(upsi, g, -urho, n, dr)
    
    erep=0.0
#!###########################################################
#!   TODO: ENTER APPROPRIATE EXPRESSION FOR
#!         REPULSION ENERGY HERE
#!###########################################################
    for i in np.arange(0,n-1,1):
        erep = erep + urho[i] * upot[i]

    print('Electron repulsion energy:',erep)
    
    for i in np.arange(0,n-1,1):
        pot[i]=upot[i]/r[i]
    
    return pot


#######################################################################################
############################ Schroedinger-solver related functions  ####################
#######################################################################################

#######################################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!
#!   Calculate function g(x) for radial Schroedinger equation in 
#!   the way it is needed for the numerow method (depends on Z, E and \phi(r) ) )
#!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#######################################################################################

def calc_g(e,z,r,pot,n):
    g=np.zeros(n)
#!###########################################################
#!   TODO: CALCULATE THE FUNCTION G(X) FOR THE RADIAL SCHROEDINGER EQUATION AS NEEDED BY THE NUMEROW METHOD
#!###########################################################
    for i in range(n):
        g[i] = 2 * (e+z/r[i]-pot[i])
    return g

############################### Initialize u(r) at core to solve Schroedinger equation #####

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!
#!  Set proper initial values of U at the boundary
#!  as needed by the Numerow method to solve the Schroedinger equation.
#!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def init_u_for_se(r,z,n):
    u=np.zeros(n)
#!###########################################################
#!   TODO: SET PROPEER INITIAL VALUES FOR THE WAVEFUNCTION (U)
#!         NEEDED BY NUMEROW METHOD (2 lines)
#!###########################################################
    for i in range(n):
        u[i] = r[i]*(1-z*r[i])
    return u
    
############################### Normalize u(r) ################################################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!
#!  Normalize U properly
#!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def normalize_u(u,r,z,n,dr):
    norm=0
#!###########################################################
#!   TODO: Complete this routine to normalize U() (approximately 3 lines)
#!###########################################################
    for i in range(n):
        norm = norm + dr*u[i]**2
    
    return u / np.sqrt(norm)

#################################### Count nodes in function u ################################

def count_nodes(u,n):
    nodes=0
#!###########################################################
#!   TODO: Complete this function to compute the number of nodes of the function U() (approximately 3 lines)
#!###########################################################
    for i in range(0, n-1):
        if (u[i] * u[i+1] < 0):
            nodes = nodes + 1
    return nodes


############################ Solve radial Schroedinger equation for a given potential ##########
############################ emin and emax specifiy energy window  ##########

def solve_se(z,pot,n,tol,r,dr):
############# initialize energy window in which the ground state energy should be searched for
    emax=0.0
    emin=-20.0
# arrays used for numerow method
    s=np.zeros(n)
    g=np.zeros(n)
# initialize radial wavefunction
    u=init_u_for_se(r,z,n)
# guess for energy of wavefunction
    e=(emax+emin)/2.0
    
#!###########################################################
#!   TODO: CONSTRUCT A WHILE LOOP THAT FINDS THE SOLUTION TO THE SCHROEDINGER EQUATION:
#
#!         CALL THE REQUIRED ROUTINES IN THE CORRECT ORDER
#!         TO OBTAIN THE SOLUTION FOR THE SCHROEDINGER EQUATION FOR EACH TRIAL ENERGY
#!         (about 6 lines)
#!
#!         ONCE THE WAVEFUNCTIOM HAS BEEN OBTAINED COMPUTE ITS NUMBER OF NODES AND
#!         MODIFY EMIN and EMAX ACCORDINGLY. COMPUTE NEW TRIAL ENERGY
#!         (about 4 lines lines)
#!###########################################################
    zs = 0
    while (np.abs(e-zs) > tol):
        zs = e
        u = init_u_for_se(r, z, n)
        g = calc_g(e, z, r, pot, n)
        u = numerow_from_inside(u, g, s, n, dr)
        u = normalize_u(u,r,z,n,dr)

        nodes = count_nodes(u, n)
        if nodes > 1:
            emin = e
        else:
            emax = e
        e=(emax+emin)/2.0
    print('Schroedinger equation energy=',e)
    return u



############# Define prarameters for radial grid
n =2000
rmax =7.0
rmin=rmax/(n+1)
dr=abs(rmax-rmin)/(n-1)

############# initialize radial grid points
r = np.linspace(rmin,rmax,n)
############# charge of helium core
z=2.0

############# initialize electronic radial potential to zero for He+
pot = np.zeros(n)

############# desired convergence parameter of ground state energy
tol=0.000001

############# solve radial Schroedinger equation for radial wavefunction
upsi=solve_se(z, pot, n, tol, r, dr)

erep=0.0
pot=solve_poisson(upsi, erep, r, n, dr)

for iteration in np.arange(0,10,1):
    print('Iteration ',iteration,' of self-consistent field cycle.')
    upsi=solve_se(z,pot,n,tol,r,dr)
    pot=solve_poisson(upsi,erep,r,n,dr)


plt.figure()
plt.plot(r,upsi/r)
plt.show()

