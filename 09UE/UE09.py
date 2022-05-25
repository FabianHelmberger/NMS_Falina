from re import T
import numpy as np
from scipy.constants import k
import matplotlib.pyplot as pl


# initialise random lattice
def rand_init(L):
    grid = (np.random.randint(2, size = (L,L))-0.5)*2
    return grid


def  energy (i,j, grid, spin_current):
    dim = np.shape(grid)[0]
    energy = 0
    

    for ind in range(len(Inn)):
        i = i + Inn[ind]
        j = j + Jinn[ind]

        neighbour_value = grid[ (i+dim)%dim ][ (j+dim)%dim ] 
        if neighbour_value != spin_current: 
            energy = energy - 1 
        else:
            energy = energy + 1 
    return J*energy

def sweep(grid, N, Temp):
    for n in range(N):
        i = round(random.random()*(L-1))
        j = round(random.random()*(L-1))
        spin_current = grid[i][j]

        E_current = energy(i, j, grid, spin_current)
        E_flipped = energy(i, j, grid, spin_current*(-1))

        
        delta_E = E_flipped - E_current
        r = np.exp(-delta_E/(k_b*Temp))

        if random.random() < min(1,r): 
            grid[i][j] = (-1) * grid[i][j]
    
    return grid


#main
L = 20              # Lattice size (i.e. NMAX x NMAX) 
J = -1.0            # Ising lattice coupling constant 
			        # J > 0: Ferromagnetic 
			        # J < 0: Antiferromagnetic
H = 0.0             # Ising lattice field strength 
Inn = np.array([0, 0, 1, -1])
Jinn = np.array([1, -1, 0, 0])
random = np.random.default_rng()
Temps = [1,2,3,4]
k_b = k
Anzahl_Gitterplätze = L*L
grid = rand_init(L)
mag_T = np.zeros( len(Temps) )
mag_square_T = np.zeros( len(Temps) )
mag_fourth_T = np.zeros( len(Temps) )
size = np.array([4,8])

for ind in range(len(size)):
    N_warm_up = size[ind]
    N_mess = size[ind]

    for temperatures in range(len(Temps)):
        Temp = Temps[temperatures]
        mag = 0 
        mag_square = 0 
        mag_fourth = 0 

        # warm up sweeps
        for n in range(N_warm_up):
            grid = sweep(grid, Anzahl_Gitterplätze, Temp)
        
        # messungen
        for n in range(N_mess):
            grid = sweep(grid, Anzahl_Gitterplätze, Temp)
            mag = mag + np.sum(grid)/(L*L)
            mag_square = mag_square + (np.sum(grid)/(L*L))**2 
            mag_fourth = mag_fourth + (np.sum(grid)/(L*L))**4

        mag_T[temperatures] = mag/N_mess
        mag_square_T[temperatures] = mag_square/N_mess
        mag_fourth_T[temperatures] = mag_fourth/N_mess

    pl.plot(mag_T)
    pl.show()


   
 

 




