from backpack import mapback as mb
from backpack import writepack as wp
import os


parent = 'lammps_protocols'
filename = 'lammps-AT-config'
path = os.path.join(parent, filename)

initial_data = wp.read_lammps_2(filename)

backmapped = mb.backmap(initial_data, 50)

coord = mb.reconfig_frame(backmapped)

hilo = [[-23.471697,23.471697],[-23.471697,23.471697],[-23.471697,23.471697]]
wp.write_lammps_input('lammps-AT-config-2',coord, hilo, 17.0)

#os.system can be used to put things on the console