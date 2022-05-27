from re import T
import numpy as np
from scipy.constants import k
import matplotlib.pyplot as pl

Dimension = 2
Next_neighbours = 2 * Dimension

Relative_neighbours = [[1, 0], [-1, 0], [0, 1], [0, -1]]

# initialise random lattice
def rand_init(L):
    grid = (np.random.randint(2, size = (L,L))-0.5)*2
    return grid


def energy (x, y, grid):
    dim_x = np.shape(grid)[0]
    dim_y = np.shape(grid)[1]
    energy = 0

    my_spin = grid[x,y]
    for rn in Relative_neighbours:
        x, y = x + rn[0], y + rn[1]
        equal = (my_spin == grid[(dim_x + x)%dim_x][(dim_y + y)%dim_y])

        if equal    : 
            energy = energy - 1 
        if not equal:
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
random = np.random.default_rng()


# Parameters than can be changed

L = 20              # Lattice size 
size = np.array([4,8])
Temps = [1,2,3,4]
J = -1.0            # Ising lattice coupling constant 
			        # J > 0: Ferromagnetic 
			        # J < 0: Antiferromagnetic


# Initialisation for code
H = 0.0             # Ising lattice field strength 
k_b = k
Anzahl_Gitterplätze = L*L
Inn = np.array([0, 0, 1, -1])
Jinn = np.array([1, -1, 0, 0])
grid = rand_init(L)
mag_T = np.zeros( len(Temps) )
mag_square_T = np.zeros( len(Temps) )
mag_fourth_T = np.zeros( len(Temps) )
binder = np.zeros( len(Temps) )


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
        # Binder Kumulante:
        binder[temperatures] = 1 - 1/3 * mag_square_T[temperatures]/(mag_square_T[temperatures] * mag_square_T[temperatures])

    pl.plot(mag_T)
    pl.title("magnetisation over temp; System size = " +str(N_mess)+   "x" + str(N_warm_up))
    pl.xlabel("Temperature")
    pl.ylabel("magnetisation")
    pl.show()

    pl.plot(mag_square_T)
    pl.title("magnetisation squared over temp; System size = " +str(N_mess)+   "x" + str(N_warm_up))
    pl.xlabel("Temperature")
    pl.ylabel("magnetisation squared")
    pl.show()

    pl.plot(mag_fourth_T)
    pl.title("magnetisation to the fourth power over temp; System size = " +str(N_mess)+   "x" + str(N_warm_up))
    pl.xlabel("Temperature")
    pl.ylabel("magnetisation^4")
    pl.show()

    pl.plot(binder)
    pl.title("Binder-Kumulante; System size = " +str(N_mess)+   "x" + str(N_warm_up))
    pl.xlabel("Temperature")
    pl.ylabel("Binder-Kumulante")
    pl.show()


   
 

 




