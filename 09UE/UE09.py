from tempfile import TemporaryFile
import numpy as np
from scipy.constants import k
import matplotlib.pyplot as pl
import copy 

k = 1
Relative_neighbours = [[1, 0], [-1, 0], [0, 1], [0, -1]]

# initialise random lattice
def rand_init(L):
    grid = (np.random.randint(2, size = (L,L))-0.5)*2
    return grid

def flip(grid_, x, y):
    grid = copy.deepcopy(grid_)
    grid[x,y] = -1*grid_[x,y]
    return grid


def energy (x, y, grid, j):
    dim = np.shape(grid)[0]
    energy = 0
    x0, y0 = x, y
    spin = grid[x0, y0]

    for coord in Relative_neighbours:
        x, y = x0 + coord[0], y0 + coord[1]
        equal_spin = (spin == grid[(dim + x)%dim][(dim + y)%dim])

        if equal_spin   : 
            energy = energy - 1 
        else            :
            energy = energy + 1 

    return j*energy


def sweep(grid, N, temp):
    sz = np.shape(grid)[0]

    for n in range(N):
        rand = np.random.randint(sz, size=2)
        x, y = rand[0], rand[1]
        E_current = energy(x, y, grid, 1              )
        E_flipped = energy(x, y, flip(grid, x, y),1   )
        delta_E = E_flipped - E_current
        r = np.exp(-delta_E/(k*temp))

        if np.random.random() < min(1,r): 
            grid = flip(grid, x, y)

    return grid

#main
random = np.random.default_rng()

# Parameters than can be changed
Sizes   = np.array([4,8,16,32])
Temps   = [i for i in range(1, 10)]
N_warm  = int(10**1)
N_meas  = int(10**2)
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

    for temp, j in zip(Temps, range(len(Temps))):
        #Temp        = Temps[j]
        mag         = 0 
        mag_square  = 0 
        mag_fourth  = 0
        print('\n')
        print("Temperature  : ", temp)
        print("Size         : ", size, 'x', size)

        # warm up sweeps
        for n in range(N_warm):
            
            print("Warmup: ",round((n+1)/(N_warm) * 100, 1),"%", end = "\r")
            grid = sweep(grid, Num, temp)
        print('\n', end = "\r")

        # measurements
        for n in range(N_meas):
            print("Measurement: ",round((n+1)/(N_meas) * 100, 1),"%", end = "\r")
            grid        = sweep(grid, Num, temp)
            mag         = mag + np.sum(grid)/Num
            mag_square  = mag_square + (np.sum(grid)/Num)**2 
            mag_fourth  = mag_fourth + (np.sum(grid)/Num)**4

        # safe grid to file
        outfile = './data/'+str(size)+'x'+str(size)+'_'+str(temp)+'.npy'

        np.save(outfile, grid)

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


   
 

 




