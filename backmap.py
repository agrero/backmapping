from backpack import mapback as mb
import pandas as pd
import os
import random

# notes
# why not have a python file that can write a lammps protocol? like duh haha

components = ['LAMMPS data file via write_data, version 7 Jan 2022, timestep = 80842034',
              'Atoms # molecular', 'Velocities', 'Bonds', 'Angles', 'Dihedrals']

initial_data = mb.read_lammps('outdata', components)

backmapped = mb.backmap(initial_data)

bonds = mb.make_bonds(backmapped)
coord = mb.reconfig_frame(backmapped)
angles = mb.make_angles(backmapped)
dihedrals = mb.make_dihedrals(backmapped)

hilo = [[-23.471697,23.471697],[-23.471697,23.471697],[-23.471697,23.471697]]
mb.write_lammps_input('lammps-AT-config',coord, bonds, angles, dihedrals, hilo, 17.0)

#os.system can be used to put things on the console