
import pandas as pd
import os 
import numpy as np

from backpack import mapback as mb
from backpack import writepack as wp
from backpack import mapanly as ma

# code here is for talapas to run the bond angle and dihedral stuff

"""chain = range(1, len(xyz[-1]) + 1)

multidex = [(i,i) for i in chain]
multidex = pd.MultiIndex.from_tuples(multidex, names=['atom','chain'])"""
"""frame = pd.DataFrame(xyz[-1], columns=['x','y','z'], index=multidex)"""
parent = 'lammps_protocols'
file = 'lammps-AT-config'

data = pd.read_csv('atomistic.csv')
data.drop(axis=1,columns=['Unnamed: 0'], inplace=True)

# 30 monomers per 500 chains
atoms= [i for i in range(1,len(data) + 1)]
chain = [i for i in range(1,len(atoms)+1) for x in range(30)]
m1 = [(atoms[i],chain[i]) for i in range(0, len(atoms))]

multi_dex = pd.MultiIndex.from_tuples(m1, names=['atom','chain'])

new_data = pd.DataFrame(data.values, index=multi_dex, columns=['x','y','z'])
path = os.path.join('lammps_protocols', 'lammps-AT-config')
coordinates = wp.read_lammps_2(path)






# uncomment these to download good plots
ma.plot_bondlength(new_data, filter_no=100, bin_width=0.05, xtick_dist=5)
ma.plot_bondlength(coordinates, filter_no=5, bin_width=0.2, xtick_dist=5)

#ma.plot_bondangle(new_data, bin_width=1, xtick_dist=8)
#ma.plot_bondangle(coordinates, bin_width=1.5, xtick_dist=10.2)

#ma.plot_dihedral(new_data, bin_width=5, xtick_dist=10)
#ma.plot_dihedral(coordinates, bin_width=5, xtick_dist=10)
