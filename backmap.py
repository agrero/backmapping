from sympy import *
from backmap import backmap as bm
import numpy as np
import math as math
import time
import pandas as pd
import random

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import timeit

def make_bonds(coordinates, no_monomers):
    #this is a monster but exactly what I want
    chains = coordinates.index.get_level_values('chain').drop_duplicates()

    bond_rep = []
    multidex_rep = []
    for i in chains:

        #chain | (atom_i, atom_j) 
        #should make a chain of 3 and 4 atoms and test to see if this method makes the correct dataframe
        bonds_list = [(i, i+1) for i in range(1, len(coordinates.loc[i]))]
        bonds = pd.DataFrame(data=bonds_list, index=[i])
        bond_rep.append(bonds)

    return pd.concat(bond_rep)

#should make a function that can take in a set of coordinates, chain number, 
rand_seed = 69420
data = {
    'x':[1.0,1.05,1.2],
    'y':[1.0,1.05,1.2],
    'z':[1.0,1.05,1.2]
}

init_coordinate = pd.DataFrame(data, columns=data.keys())
#init_coordinate.index.name = 'chain'

radius = 1
volume = np.pi * radius ** 2




"""
df2 = pd.DataFrame(init_coordinate)
print(df2)
fudged = pd.DataFrame(bm.fudge_position(df2, radius=radius))
print('fudged')
print(fudged)

multi_dex_1, multi_dex_2 = bm.gen_multidex(range(len(df2)))

merged = pd.concat([
    pd.DataFrame(df2.values, multi_dex_1),
    pd.DataFrame(fudged.values, multi_dex_2)
])
merged.sort_index(level='chain', inplace=True)
merged.rename(columns={0:'x', 1:'y', 2:'z'}, inplace=True)
print('merged \n')
print(merged)
fudged = pd.DataFrame(bm.fudge_position(merged, radius=1))
print('fudged 2\n')
print(fudged)
"""

"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(merged.iloc[:,0], merged.iloc[:,1], merged.iloc[:,2], s=volume)
ax.view_init(30,185)
plt.show()
"""
