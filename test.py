
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

data = np.load('ALLXYZ.npy')
frame = pd.DataFrame(data[-1][0])
frame.to_csv('atomistic.csv')
values = []
for i in data[-1]:
    frame = pd.DataFrame(i)
    values.append(list(i))

value_f = sum(values, [])
value_frame = pd.DataFrame(value_f)
value_frame.to_csv('atomistic.csv')
#writename = os.path.join(parent, file)
#coordinates = wp.read_lammps_2(writename)

#ma.plot_bondlength(coordinates, filter_no=2, bin_width=.1)
#ma.plot_bondangle(coordinates, bin_width=5)
#ma.plot_dihedral(coordinates, bin_width=5)
