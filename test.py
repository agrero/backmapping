from backpack import mapback as mb
from backpack import writepack as wp
import default_configs as dc

import os 

"""import pandas as pd

initial_data = {
    'x': [1.0001, 2.0001, 3.0001],
    'y': [1.0001, 2.0001, 3.0001],
    'z': [1.0001, 2.0001, 3.0001]
}
initial_data = pd.DataFrame(initial_data)

backmapped = mb.backmap(initial_data)
for i in range(6):
    backmapped = mb.backmap(backmapped)
    """
file_1 = 'lammps-AT-config'
file_2 = 'outdata'

initial_input = os.path.join('lammps_protocols', file_2)

pastor = wp.read_lammps_2(initial_input)

print(pastor)