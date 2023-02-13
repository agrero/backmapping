from sympy import *

import numpy as np
import math as math
import time
import pandas as pd
import random

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def fudge_position(data, radius):
    """
    the radius being divided by 100 will need to change
    """
    #fudgers = [[random.uniform(-radius/500, radius/500) for i in range(3)] for i in range(len(data))]
    fudgers = [[radius/10 for i in range(3)] for i in range(len(data))]
    fudgers = pd.DataFrame(fudgers, columns=['x','y','z'])

    new_data = data + fudgers

    return new_data

def gen_multidex(no_dex):

    #i'm just going to duplicate the code below but should try to be more clever later
    m1 = [(x,1) for x in range(len(no_dex))]
    m2 = [(x,2) for x in range(len(no_dex))]
    
    multi_dex_1 = pd.MultiIndex.from_tuples(m1, names=['chain', 'atom'])
    multi_dex_2 = pd.MultiIndex.from_tuples(m2, names=['chain', 'atom'])

    return multi_dex_1, multi_dex_2

def make_bonds(coordinates, no_monomers):
    #this is a monster but exactly what I want
    chains = coordinates.index.get_level_values('chain').drop_duplicates()

    for i in chains:
        query
    
    


rand_seed = 69420
data = {
    'x':[1.0,1.05,1.2],
    'y':[1.0,1.05,1.2],
    'z':[1.0,1.05,1.2]
}

init_coordinate = pd.DataFrame(data, columns=data.keys())
radius = 5
volume = np.pi * radius ** 2


multi_dex_1, multi_dex_2 = gen_multidex(range(len(init_coordinate)))
df2 = pd.DataFrame(init_coordinate)
fudged = pd.DataFrame(fudge_position(df2, radius=radius))

multi_dex_1, multi_dex_2 = gen_multidex(range(len(df2)))

merged = pd.concat([
    pd.DataFrame(df2.values, multi_dex_1),
    pd.DataFrame(fudged.values, multi_dex_2)
])
merged.sort_index(level='chain', inplace=True)

a = make_bonds(merged, 2)


"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(merged.iloc[:,0], merged.iloc[:,1], merged.iloc[:,2], s=volume)
ax.view_init(30,185)
plt.show()
"""