def write_lammps_config(config_list, filename):
    #can probably rewrite this as a dictionary
    with open(filename, 'w') as f:
        for i in config_list:
            f.write(f'{i}\n')

#######
# Parameters
#######
# should likely distinguish between if you're reading a restart 
# or you're reading an input file
# need a way to change the mass in here CANNNOT forget to do that 
param_dict ={
    'atom_style' : 'molecular', 
    'boundary': 'p p p',
    'units' : 'real',
    'read_data' : 'lammps-AT-config',
    'read_restart' : 'RESTART',
    'bond_style' : 'harmonic',
    'bond_coeff' : '1 450 1.54',
    'special_bonds' : 'lj 0.0 0.0 0.0 angle no dihedral no',
    'angle_style' : 'harmonic',
    'angle_coeff' : '1 61.875 114.00',
    'dihedral_style' : 'opls',
    'dihedral_coeff' : '* 1.14110, -0.27084 3.143 0.0',
    'mass' : '* 14.0',
    'pair_style' : 'lj/mdf 12.0 14.0',
    'pair_coeff' : '1 1 0 0.0912 2.95 12.0 14.0',
    'neigh_modify' : 'every 1 delay 0 cehck yes',
    'thermo_style' : 'custom step temp press ke pe etotal',
    'thermo_modify' : 'lost warn',
    'dump' : 'coords all atom 1 dump.xyz.dat',
    'dump_modify' : 'coords',
    'restart' : '10000 restart.dat'
} 
# NOT INCLUDED
    # minimize
    # run
    # timestep
# this is how it should write the parameters
# we should also not forget to ask if you'd want to write the data at each step
# in order to make a video!
for key in param_dict:
    print(f'{key} {param_dict[key]}')

atom_style = 'atom_style molecular'
boundary = 'boundary p p p'
units = 'units real'

read_data = 'read_data lammps-AT-config'
read_restart = 'read_restart RESTART'

bond_style = 'bond_style harmonic'
bond_coeff = 'bond_coeff 1 450 1.54'
special_bonds = 'special_bonds lj 0.0 0.0 0.0 angle no dihedral no'

angle_style = 'angle_style harmonic'
angle_coeff = 'angle_coeff 1 61.875 114.00'

dihedral_style = 'dihedral_style opls'
dihedral_coeff = 'dihedral_coeff * 1.14110 -0.27084 3.143 0.0'

mass = 'mass * 14.0'

pair_style = 'pair_style lj/mdf 12.0 14.0'
pair_coeff = 'pair_coeff 1 1 0.0912 3.95 12.0 14.0'

neigh_modify = 'neigh_modify every 1 delay 0 check yes'
minimize = 'minimize 1.0e-4 1.0e-6 1000 1000'

thermo_style = 'thermo_style custom step temp press ke pe etotal'
thermo_modify = 'thermo_modify lost warn'

dump = 'dump coords all atom 1 dump.xyz.dat'
dump_modify = 'dump_modify coords scale no'

restart = 'restart 10000 restart.dat'

#######
# Run Conditions
#######

timestep = 'timestep 0.02'
run = 'run 100'
#ask if you wanna minimize between each run
# minimize is the same for each
# i think the beginning parameters can be put in a dict but the run 
# should be written through a function

timestep_2 = 'timestep 0.75'
run_2 = 'run 100'
#minimize

timestep_3 = 'timestep 1.0'
run_3 = 'run 100'
#minimize

timestep_4 = 'timestep 1.25'
run_4 = 'run 20000'

write_restart = 'write_restart FINALCONFIG'
stuff_list = [atom_style, boundary, units, read_restart, bond_style, bond_coeff, special_bonds,
              angle_style, angle_coeff, dihedral_style, dihedral_coeff, mass, pair_style,
              pair_coeff, neigh_modify, thermo_style, thermo_modify, dump, dump_modify, restart, timestep,run, minimize,
              timestep_2, run_2, minimize,timestep_3, run_3, minimize,timestep_4, run_4, write_restart]

#write_lammps_config(stuff_list, 'lammps_equil_config')