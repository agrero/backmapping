
import pandas as pd
import os 
import numpy as np

from backpack import mapback as mb
from backpack import writepack as wp


# code here is for talapas to run the bond angle and dihedral stuff

xyz = wp.read_lammps_2('lammps-AT-config')
frame = pd.DataFrame(xyz)
print(frame)
"""chain = range(1, len(xyz[-1]) + 1)

multidex = [(i,i) for i in chain]
multidex = pd.MultiIndex.from_tuples(multidex, names=['atom','chain'])"""
"""frame = pd.DataFrame(xyz[-1], columns=['x','y','z'], index=multidex)"""



hilo=[[-50,50],[-50,50],[-50,50]]
# max in the frame is 50
wp.write_lammps_input('allcgxyz.lmp', coordinates=frame, hi_lo=hilo, mass=112.0)