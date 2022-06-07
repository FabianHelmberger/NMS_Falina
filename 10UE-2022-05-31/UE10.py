from sys import stderr
import numpy as np
from sympy import Lambda
import copy
import matplotlib.pyplot as plt
import logging
logging.basicConfig(filename="./data/log.txt")

# read data
Raw_data = np.loadtxt("SimMRimage.dat", dtype=int)
X_max = (Raw_data[-1, 0])
Y_max = (Raw_data[-1, 1])
Data = np.reshape(Raw_data[:, 2], (X_max, Y_max))

# Constants
Lambda_ = 1.5
T_f = 0.1
T_i = 10
J = 1.5
Num = X_max * Y_max
N_wamrup = 2
Relative_neighbours = [[1, 0], [-1, 0], [0, 1], [0, -1]]

# stats
Stats = {   '1'    : [30   , 30    , 'BG' ], 
            '2'    : [426  , 59    , 'WM' ],
            '3'    : [602  , 102   , 'GM' ],
            '4'    : [1223 , 307   , 'CSF'],
            '5'    : [167  , 69    , 'SB' ]}


def mutate(grid_, x, y):
    grid = copy.deepcopy(grid_)
    new = np.random.randint(1,6)

    while (new == grid[x,y]):
        new = np.random.randint(1,6)
        
    grid[x,y] = new
    return grid

def energy(nmr, Labels, pos):
    energy = 0
    x0, y0 = pos[0], pos[1]

    #for x0, y0 in zip(range(X_max), range(Y_max)):
    label  = Labels[x0, y0]

    for coord in Relative_neighbours:
        x, y = x0 + coord[0], y0 + coord[1]
        equal_label = (label == Labels[(X_max + x)%X_max][(Y_max + y)%Y_max])

        if equal_label   :
            energy = energy - J
        
    mean = Stats[str(label)][0]
    std  = Stats[str(label)][1]
    int = nmr[x0, y0]
    energy = energy + (int - mean)**2/(2*std**2) + np.log(std)

    return energy

def sweep(nmr, Labels, T):

    grid = Labels
    for n in range(Num):
        
        rand = np.random.randint(1, [X_max, Y_max])
        x, y = rand[0], rand[1]
        grid_mut    = mutate(grid, x, y)
        E_current   = energy(nmr, grid, [x,y]   )
        E_mutated   = energy(nmr, grid_mut, [x,y])
        delta_E     = E_mutated - E_current

        if np.random.random() < min(1,np.exp(-delta_E/T)): 
            grid = grid_mut
        print("sweep: ",round((n+1)/(Num) * 100, 1),"%", end = "\r")
    print('\n')
    return grid


# Main
# initialise a random integer grid
grid = np.random.randint(1, 6, size=(X_max, Y_max))

temps   = [T_i]
T       = T_i

# Warm-up loop
for i in range(N_wamrup):
    grid = sweep(Data, grid, T)

while (T > T_f):
    T = T / Lambda_

    print('now starting T = '+str(T))
    grid = sweep(Data, grid, T)
    outfile = './data/'+str(T)+'.npy'
    np.save(outfile, grid)
    temps.append(T)

logging.debug('J = ',J)
logging.debug('Temps = ')

for T in temps:
    logging.debug(T)

plt.pcolormesh(grid)
plt.show()