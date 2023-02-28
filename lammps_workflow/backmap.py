import talapas_backmapping as bm
import pandas as pd
import os

components = ['LAMMPS data file via write_data, version 7 Jan 2022, timestep = 80842034',
              'Atoms # molecular', 'Velocities', 'Bonds', 'Angles', 'Dihedrals']

os.chdir(os.path.join(os.getcwd(), 'lammps_workflow'))

initial_data = bm.read_lammps('outdata.txt', components)
backmapped = bm.backmap(initial_data)

bonds = bm.make_bonds(backmapped)
coord = bm.reconfig_frame(coordinates=backmapped)
angles = bm.make_angles(backmapped)
dihedrals = bm.make_dihedrals(backmapped)


hilo = [[-23.471697,23.471697],[-23.471697,23.471697],[-23.471697,23.471697]]
#bm.write_lammps_input('test.txt',coord, bonds, angles, dihedrals, hilo, 17.0)