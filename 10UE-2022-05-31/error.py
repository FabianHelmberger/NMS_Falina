import numpy as np

#data = np.load('./6.666666666666667.npy')
#print(data)

import os
data = []
for filename in os.listdir(os.getcwd()+"/data"):
    temp = np.load(os.getcwd()+"/data/"+filename)
    data.append(temp)
print(np.shape(np.array(data)))