import numpy as np
import matplotlib.pyplot as plt
from sympy import E
import math

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
 
def numerow_from_inside(u,g,s,n,dr):
#!###########################################################
#!   TODO: ENTER APPROPRIATE NUMEROW EXPRESSION HERE
#!###########################################################
    for i in range(1,n-1):
        u[i+1] = (1/(1+((dr**2)/12)*g[i+1]))*((2*u[i]*(1-((5/12)*(dr**2)*g[i])))-(u[i-1]*(1+(((dr**2)/12)*g[i-1])))+(((dr**2)/12)*(s[i+1]+10*s[i]+s[i-1])))   
    return u   

def numerow_from_outside(u,g,s,n,dr):
#!###########################################################
#!   TODO: ENTER APPROPRIATE NUMEROW EXPRESSION HERE
#!###########################################################
    for i in range(n-2,0,-1):
       u[i-1] = (1/(1+((dr**2)/12)*g[i-1]))*((2*u[i]*(1-((5/12)*(dr**2)*g[i])))-(u[i+1]*(1+(((dr**2)/12)*g[i+1])))+(((dr**2)/12)*(s[i+1]+10*s[i]+s[i-1])))   
    return u   


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
        urho[i] = u[i]**2
    return urho


def solve_poisson(upsi,erep,r,n,dr):
# arrays used for numerow method

    # array needed for numerow method
    g=np.zeros(n)
    s=np.zeros(n)
    upot=np.zeros(n)
#!###########################################################
#!   TODO: CALL THE THREE REQUIRED FUNCTIONS IN THE CORRECT ORDER
#!         TO OBTAIN THE SOLUTION FOR THE POISSON EQUATION HERE (3 lines)
#!###########################################################
    urho=calc_urho(upsi,r,n)
    
    for i in range(n):
        s[i]=-4.*math.pi*urho[i]*r[i]
        
    upot=init_upot(r,n)
    upot=numerow_from_outside(upot,g,s,n,dr)
    for i in range(n):
        upot[i]=upot[i]/r[i]

    erep=0.0
#!###########################################################
#!   TODO: ENTER APPROPRIATE EXPRESSION FOR
#!         REPULSION ENERGY HERE
#!###########################################################
    f=np.zeros(n)
    for i in range(0,n):
        f[i]=4.*math.pi*r[i]*r[i]*urho[i]*upot[i]
               
    erep=trapezoid(f,r,dr)
    
    print('Electron repulsion energy:',erep)
    
    return upot,erep


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
from scipy.integrate import trapezoid
from numpy.linalg import norm

def normalize_u(u):
#!###########################################################
#!   TODO: Complete this routine to normalize U() (approximately 3 lines)
#!###########################################################
    f=np.zeros(n)
    for i in range(n):
        f[i]=4*math.pi*(r[i]*r[i])*u[i]*u[i]
    
    norm=trapezoid(f,r,dr)
    return u/np.sqrt(norm)


#################################### Count nodes in function u ################################

def count_nodes(u,n):
    nodes=0
#!###########################################################
#!   TODO: Complete this function to compute the number of nodes of the function U() (approximately 3 lines)
#!###########################################################
    for i in range(n-1):
        if (u[i] * u[i+1] < 0):
            nodes = nodes + 1
    return nodes


############################ Solve radial Schroedinger equation for a given potential ##########
############################ emin and emax specifiy energy window  ##########

def solve_se(z,pot,n,tol, r):
############# initialize energy window in which the ground state energy should be searched for
    emax=0.0
    emin=-20.0
    e=(emax+emin)/2
# arrays used for numerow method
    s=np.zeros(n)
    g=calc_g(e,z,r,pot,n)
#how many iterations where neccessary
    Iterations=0
    
# initialize radial wavefunction

    u_init=init_u_for_se(r,z,n)
    #print(dir(u))
    u=numerow_from_inside(u_init,g,s,n,dr)
    
    nodes=count_nodes(u,n)
    
    for i in range(1,100):
        if count_nodes(u,n)>0:
            emax = e
            e = (emax + emin) / 2.0
            g = calc_g(e, z, r, pot, n)
            u = numerow_from_inside(u, g, s, n, dr)
            #print(count_nodes(u, n), emin)
        elif count_nodes(u,n)==0:
            emin = e
            e = (emax + emin) / 2.0
            g = calc_g(e, z, r, pot, n)
            u = numerow_from_inside(u, g, s, n, dr)
            #print(count_nodes(u, n), emin
          
        #print(nodes)
           
        Iterations=Iterations+1   
        
    psi=u/r    
    emax=e
    e=(emax+emin)/2.0
    
    print('Groundstate Energy=',e)
    
    psi_norm=normalize_u(psi)
    
    return psi_norm,e



############# Define prarameters for radial grid
n =2000
rmax =7.0
rmin=rmax/(n+1)
dr=abs(rmax-rmin)/(n-1)

############# initialize radial grid points
r = np.linspace(rmin,rmax,n)

############# charge of helium core
z=2.0




# Question 1
############# initialize electronic radial potential to zero for He+
pot = np.zeros(n)

############# desired convergence parameter of ground state energy
tol=0.000001

############# solve radial Schroedinger equation for radial wavefunction
print('now solving SE for He...')
upsi    =   solve_se(z,pot,n,tol,r)[0]
e       =   solve_se(z,pot,n,tol,r)[1]

#plt.plot(r,upsi,label='$\Psi$')
#plt.title('Upsi for Pot=0')
#plt.legend()
#plt.show()


# Question 2
print('now solving PE for HE...')
erep    =   0.0
upot    =   solve_poisson(upsi, erep, r, n, dr)[0]
erep    =   solve_poisson(upsi, erep, r, n, dr)[1]

#plt.plot(r,upot,label='$\Phi$')
#plt.plot(r,upsi,label='$\Psi$')
#plt.legend()
#plt.show()

# Question 3
upsi    =   solve_se(z,pot,n,tol,r)[0]
e       =   solve_se(z,pot,n,tol,r)[1]
e_1     =   np.inf


pot     =   solve_poisson(upsi,erep,r,n,dr)[0]
erep    =   solve_poisson(upsi,erep,r,n,dr)[1]
iter    =   0
print("\n")
print("starting loop for HE")
while(np.abs(e-e_1)>tol):

    e_1 =   e
    print('Iteration ',iter,' of self-consistent field cycle.')
    upsi    =   solve_se(z,pot,n,tol,r)[0]
    e       =   solve_se(z,pot,n,tol,r)[1]
    pot     =   solve_poisson(upsi,erep,r,n,dr)[0]
    erep    =   solve_poisson(upsi,erep,r,n,dr)[1]
    iter=iter+1

E_ges   =   2*e-erep

print('Groundstate He:',E_ges,'E_H')

upot    =   solve_poisson(upsi,erep,r,n,dr)[0]
erep    =   solve_poisson(upsi,erep,r,n,dr)[1]

#plt.plot(r,upot,label='$\Phi$')
#plt.plot(r,upsi,label='$\Psi$')
#plt.legend()
#plt.plot(r,upsi/r)
#plt.show()




