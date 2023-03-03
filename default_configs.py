from backpack import writepack as wp

#######
# Parameters for writing lammps config
#######
default_peequil_dict ={
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
    'neigh_modify' : 'every 1 delay 0 check yes',
    'minimize' : '1.0e-4 1.0e-6 1000 1000',
    'thermo_style' : 'custom step temp press ke pe etotal',
    'thermo_modify' : 'lost warn',
    'dump' : 'coords all atom 1 dump.xyz.dat',
    'dump_modify' : 'coords',
    'restart' : '10000 restart.dat'
} 

default_soft_dict = {
    'atom_style' : 'molecular',
    'boundary' : 'p p p',
    'units' : 'real',
    'read_data' : 'lammps-AT-config',
    'bond_style' : 'harmonic',
    'bond_coeff' : '1 450 1.54',
    'special_bonds' : 'lj 0.0 0.0 0.0 angle no dihedral no',
    'angle_style' : 'harmonic',
    'angle_coeff' : '1 450 1.54',
    'dihedral_style' : 'opls',
    'dihedral_coeff' : '* 1.14110 -0.27084 3.143 0.0',
    'mass' : '* 14.0',
    'pair_style' : 'soft 10.5',
    'pair_coeff' : '* * 10.0',
    'variable' : 'prefactor equal ramp(10,100)',
    'fix' : '1 all adapt 1 pair soft a * * v_prefactor',
    'neighbor' : '0.4 bin',
    'neigh_modify' : 'every 1 delay 0 check yes',
    'thermo_style' : 'custom step temp press ke pe etotal',
    'thermo_modify' : 'lost warn',
    'thermo' : '1000',
    'dump' : 'coords all atom 1 dump.xyz.dat',
    'dump_modify' : 'coords scale no',
}
# maybe later we try to do the run nos and timesteps like we did the beginning and end things 
# in the backmapping protocol writing

# we should also not forget to ask if you'd want to write the data at each step
# in order to make a video!

#######
# Inputs for writing workflow bash
#######

default_workflow_dict = {
    'job-name' : '"N10MDFrho0315"', # the second parentheses are important! 
    'output' : '"N10.%j.%N.out"',
    'partition' : 'short',
    'time' : '0-01:00:00',
    'error': 'joberr',
    'nodes' : '1',
    'ntasks-per-node' : '5',
    'cpus-per-task' : '1',
    'export' : 'ALL',
    'mail-user' : 'email@address.lll',
    'mail-type' : ['begin', 'end'],
    'account' : 'guenzagrp'
}
