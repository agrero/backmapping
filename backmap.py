from backpack import mapback as mb
from backpack import writepack as wp
import pandas as pd
import os
import random

# notes
# why not have a python file that can write a lammps protocol? like duh haha

components = ['Atoms # molecular', 'Velocities', 'Bonds', 'Angles', 'Dihedrals']

parent = 'lammps_protocols'
initial_data = wp.read_lammps('outdata', components)

backmapped = mb.backmap(initial_data)

bonds = mb.make_bonds(backmapped)
coord = mb.reconfig_frame(backmapped)
angles = mb.make_angles(backmapped)
dihedrals = mb.make_dihedrals(backmapped)

hilo = [[-23.471697,23.471697],[-23.471697,23.471697],[-23.471697,23.471697]]
mb.write_lammps_input('lammps-AT-config',coord, bonds, angles, dihedrals, hilo, 17.0)

#os.system can be used to put things on the console