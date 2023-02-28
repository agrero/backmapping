from backmap import backmapping as bm
import pandas as pd

initial_data = {
    'x': [1.0001, 2.0001, 3.0001],
    'y': [1.0001, 2.0001, 3.0001],
    'z': [1.0001, 2.0001, 3.0001]
}
initial_data = pd.DataFrame(initial_data)

backmapped = bm.backmap(initial_data)
for i in range(6):
    backmapped = bm.backmap(backmapped)
    
#print(backmapped)
bonds = bm.make_bonds(backmapped)
coord = bm.reconfig_frame(coordinates=backmapped)
angles = bm.make_angles(backmapped)
dihedrals = bm.make_dihedrals(backmapped)


hilo = [[-10,10],[-10,10],[-10,10]]
bm.write_lammps_input('test.txt',coord, bonds, angles, dihedrals, hilo, 17.0)

components = ['Atoms', 'Bonds', 'Angles', 'Dihedrals']
yes = bm.read_lammps('test.txt', components)

back_yes = bm.backmap(yes)
bonds_2 = bm.make_bonds(back_yes)
coord_2 = bm.reconfig_frame(back_yes)
angles_2 = bm.make_angles(back_yes)
dihedrals_2 = bm.make_dihedrals(back_yes)

bm.write_lammps_input('lammps_input', coord_2, bonds_2, angles_2, dihedrals_2, hilo, 17.0)
