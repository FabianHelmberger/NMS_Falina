from re import T
import numpy as np
from scipy.constants import k
import matplotlib.pyplot as pl
k = 1

Relative_neighbours = [[1, 0], [-1, 0], [0, 1], [0, -1]]

# initialise random lattice
def rand_init(L):
    grid = (np.random.randint(2, size = (L,L))-0.5)*2
    return grid

def flip(grid, x, y):
    grid[x][y] = -grid[x][y]
    return grid

def energy (x, y, grid, j):
    dim_x = np.shape(grid)[0]
    dim_y = np.shape(grid)[1]
    energy = 0
    x0, y0 = x, y
    spin = grid[x0, y0]

    for coord in Relative_neighbours:
        x, y = x0 + coord[0], y0 + coord[1]
        equal_spin = (spin == grid[(dim_x + x)%dim_x][(dim_y + y)%dim_y])
        
        #print(equal_sgrid[(dim_x + x)%dim_x][pin)
        if equal_spin   : 
            energy = energy - 1 

        if not equal_spin:
            energy = energy + 1 

    return j*energy

def sweep(grid, N, temp):
    size = np.shape(grid)[0]
    for n in range(N):
        rand = np.random.randint(size, size=2)
        x, y = rand[0], rand[1]
        print(x,y)
        #spin = grid[x][y]

        E_current = energy(x, y, grid, -1              )
        E_flipped = energy(x, y, flip(grid, x, y),-1   )

        
        delta_E = E_flipped - E_current
        r = np.exp(-delta_E/(k*temp))
        #print(r)
        if np.random.random() < np.minimum(1,r): 
            #print('flip!')
            # grid[x][y] = (-1) * grid[x][y]
            grid = flip(grid, x, y)
        else:
            print("no flip!")
    return grid


#main
random = np.random.default_rng()

# Parameters than can be changed
Sizes = np.array([4,8])
Temps = [i for i in range(1,15)]
J = -1.0            # Ising lattice coupling constant 
			        # J > 0: Ferromagnetic 
			        # J < 0: Antiferromagnetic


# Initialisation for code
H = 0.0             # Ising lattice field strength 

mag_T        = np.zeros( len(Temps) )
mag_square_T = np.zeros( len(Temps) )
mag_fourth_T = np.zeros( len(Temps) )
binder       = np.zeros( len(Temps) )

for size in Sizes:
    Num  = size**2
    grid = rand_init(size)  

    for j in range(len(Temps)):
        Temp        = Temps[j]
        mag         = 0 
        mag_square  = 0 
        mag_fourth  = 0 

        # warm up sweeps
        for n in range(size):
            grid = sweep(grid, Num, Temp)
        
        # messungen
        for n in range(size):
            grid        = sweep(grid, Num, Temp)
            #print(grid)
            mag         = mag + np.sum(grid)/Num
            mag_square  = mag_square + (np.sum(grid)/Num)**2 
            mag_fourth  = mag_fourth + (np.sum(grid)/Num)**4

        mag_T[j]        = mag/size
        mag_square_T[j] = mag_square/size
        mag_fourth_T[j] = mag_fourth/size
        # Binder Kumulante:
        binder[j] = 1 - 1/3 * mag_square_T[j]/(mag_square_T[j] * mag_square_T[j])

    pl.plot(mag_T)
    pl.title("magnetisation over temp; System size = " +str(size)+   "x" + str(size))
    pl.xlabel("Temperature")
    pl.ylabel("magnetisation")
    pl.show()

    pl.plot(mag_square_T)
    pl.title("magnetisation squared over temp; System size = " +str(size)+   "x" + str(size))
    pl.xlabel("Temperature")
    pl.ylabel("magnetisation squared")
    pl.show()

    pl.plot(mag_fourth_T)
    pl.title("magnetisation to the fourth power over temp; System size = " +str(size)+   "x" + str(size))
    pl.xlabel("Temperature")
    pl.ylabel("magnetisation^4")
    pl.show()

    pl.plot(binder)
    pl.title("Binder-Kumulante; System size = " +str(size)+   "x" + str(size))
    pl.xlabel("Temperature")
    pl.ylabel("Binder-Kumulante")
    pl.show()


   
 

 




