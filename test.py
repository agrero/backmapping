from backmap import backmapping as bm
import pandas as pd

initial_data = {
    'x': [1.0001, 2.0001, 3.0001],
    'y': [1.0001, 2.0001, 3.0001],
    'z': [1.0001, 2.0001, 3.0001]
}
initial_data = pd.DataFrame(initial_data)

backmapped = bm.backmap(initial_data)
for i in range(2):
    backmapped = bm.backmap(backmapped)
    
#print(backmapped)
bonds = bm.make_bonds(backmapped)
coord = bm.reconfig_frame(coordinates=backmapped)

angles = bm.make_angles(backmapped)
print(angles)


#hilo = [[-10,10],[-10,10],[-10,10]]
#bm.write_lammps_input(coord, hilo, 17.0)

