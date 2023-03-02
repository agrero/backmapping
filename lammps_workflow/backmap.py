import talapas_backmapping as bm
import pandas as pd
import os
import random

#notes
# why not have a python file that can write a lammps protocol? like duh haha

components = ['LAMMPS data file via write_data, version 7 Jan 2022, timestep = 80842034',
              'Atoms # molecular', 'Velocities', 'Bonds', 'Angles', 'Dihedrals']

os.chdir(os.path.join(os.getcwd(), 'lammps_workflow'))

initial_data = bm.read_lammps('outdata', components)

backmapped = bm.backmap(initial_data)

bonds = bm.make_bonds(backmapped)
coord = bm.reconfig_frame(coordinates=backmapped)

angles = bm.make_angles(backmapped)

dihedrals = bm.make_dihedrals(backmapped)



radius = 23.471697
fudgers = [[random.uniform(-radius, radius) for i in range(3)] for i in range(300)]
chain = list(range(300))
multi_dex_list = [(chain[i]+1, chain[i]+1) for i in chain]
multi_dex = pd.MultiIndex.from_tuples(multi_dex_list, names=['atom', 'chain'])
backmap = pd.DataFrame(fudgers, columns=['x','y','z'], index=multi_dex)

backmapped = bm.backmap(backmap)
coord = bm.reconfig_frame(backmapped)

bonds = bm.make_bonds(backmapped)
angles = bm.make_angles(backmapped)
dihedrals = bm.make_dihedrals(backmapped)

hilo = [[-23.471697,23.471697],[-23.471697,23.471697],[-23.471697,23.471697]]
bm.write_lammps_input('lammps-AT-config',coord, bonds, angles, dihedrals, hilo, 17.0)