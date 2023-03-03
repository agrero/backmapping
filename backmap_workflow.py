from backpack import mapback as mb
from backpack import writepack as wp
import default_configs as dc

import os

# checking for dependent directories
dir_list = os.listdir()


if not 'lammps_protocols' in dir_list:
    os.mkdir('lammps_protocols')
if not 'lammps-AT-config' in os.listdir('lammps_protocols'):
    print('no data detected, let me fix that for ya')
    test_data = mb.generate_test_data(no_monomers=30, bounds=23.471697)
    bonds = mb.make_bonds(test_data)
    angles = mb.make_angles(test_data)
    dihedrals = mb.make_dihedrals(test_data)
    mass = 17.0
    hi_lo = [[-23.4711697, 23.4711697] for i in range(3)]
    wp.write_lammps_input('lammps-AT-config', mb.reconfig_frame(test_data), 
                          bonds, angles, dihedrals, hi_lo, mass)
    
#writing protocols/configurations
bash_filename = 'bashmap.sh'
test_modules = ['python']
wp.write_backmapping_protocol(filename=bash_filename, modules=None, threads=5, no_iter=2,
                              lammps_call='start')
soft_timesteps = [1.00]
soft_run_nos = [20000]
wp.write_lammps_config(config_dict=dc.default_soft_dict, filename='in.PEsoft', 
                       restart=False, timesteps=soft_timesteps, run_nos=soft_run_nos)
equil_timesteps = [0.025,0.08,0.25,0.5,0.75,1.0,1.25]
equil_run_nos = [100,100,100,100,100,100,20000]
wp.write_lammps_config(config_dict=dc.default_peequil_dict, filename='in.PEequil', 
                       restart=True, timesteps=equil_timesteps, run_nos=equil_run_nos)

# run the protocol
os.chdir('lammps_protocols')
os.system('chmod u+x bashmap.sh')
os.system('bash bashmap.sh')