from backpack.mapback import backmapping as mb
import pandas as pd

initial_data = {
    'x': [1.0001, 2.0001, 3.0001],
    'y': [1.0001, 2.0001, 3.0001],
    'z': [1.0001, 2.0001, 3.0001]
}
initial_data = pd.DataFrame(initial_data)

backmapped = mb.backmap(initial_data)
for i in range(6):
    backmapped = mb.backmap(backmapped)
    
#print(backmapped)
bonds = mb.make_bonds(backmapped)
coord = mb.reconfig_frame(coordinates=backmapped)
angles = mb.make_angles(backmapped)
dihedrals = mb.make_dihedrals(backmapped)


hilo = [[-10,10],[-10,10],[-10,10]]
mb.write_lammps_input('test.txt',coord, bonds, angles, dihedrals, hilo, 17.0)

components = ['Atoms', 'Bonds', 'Angles', 'Dihedrals']
yes = mb.read_lammps('test.txt', components)

back_yes = mb.backmap(yes)
bonds_2 = mb.make_bonds(back_yes)
coord_2 = mb.reconfig_frame(back_yes)
angles_2 = mb.make_angles(back_yes)
dihedrals_2 = mb.make_dihedrals(back_yes)
# likely need to chagne the input here for the reading file as the 
# chain number is 384 and only one monomer and this should be changed
mb.write_lammps_input('lammps_input.txt', coord_2, bonds_2, angles_2, dihedrals_2, hilo, 17.0)
