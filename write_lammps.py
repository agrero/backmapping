def write_lammps_config(config_list, filename):
    with open(filename, 'w') as f:
        for i in config_list:
            f.write(f'{i}\n')

atom_style = 'atom_style molecular'
boundary = 'boundary p p p'
units = 'units real'

read_data = 'read_data lammps-AT-config'

bond_style = 'bond_style harmonic'
bond_coeff = 'bond_coeff 1 450 1.54'
special_bonds = 'special_bonds lj 0.0 0.0 0.0 angle no dihedral no'

angle_style = 'harmonic'
angle_coeff = 'angle_coeff 1 61.875 114.00'

dihedral_style = 'dihedral_style opls'
dihedral_coeff = 'dihedral_coeff * 1.14110 -0.27084 3.143 0.0'

mass = 'mass * 14.0'

pair_style = 'pair_style soft 10.5'
pair_coeff = 'pair_coeff * * 10.0'
variable_prefactor = 'variable prefactor equal ramp(10,100)'
fix_1 = 'fix 1 all adapt 1 pair soft a * * v_prefactor'

neighbor = 'neighbor 0.4 bin'
neigh_modify = 'neigh_modify every 1 delay 0 check yes'

thermo_style = 'thermo_style custom step temp press ke pe etotal'
thermo_modify = 'thermo_modify lost warn'
thermo = 'thermo 1000'

dump = 'dump coords all atom 1 dump.xyz.dat'
dump_modify = 'dump_modify coords scale no'

restart = 'restart 10000 restart.dat'

timestep = 'timestep 1.00'
run = 'run 20000'

write_restart = 'write_restart FINALCONFIG'
stuff_list = [atom_style, boundary, units, read_data, bond_style, bond_coeff, special_bonds,
              angle_style, angle_coeff, dihedral_style, dihedral_coeff, mass, pair_style,
              pair_coeff, variable_prefactor, fix_1, neighbor, neigh_modify, thermo_style,
              thermo_modify, thermo, dump, dump_modify, restart, run, write_restart]

write_lammps_config(stuff_list, 'lammps_config')
"""
atom_style molecular

boundary p p p

units real

read_data lammps-AT-config

bond_style harmonic
bond_coeff 1 450 1.54

special_bonds lj 0.0 0.0 0.0 angle no dihedral no

angle_style harmonic
angle_coeff 1 61.875 114.00

dihedral_style opls
dihedral_coeff * 1.4110 -0.27084 3.143 0.0

mass * 14.0

pair_style soft 10.5
pair_coeff * * 10.0
variable prefactor equal ramp(10,100)
fix 1 all adapt 1 pair soft a * * v_prefactor

neighbor        0.4 bin
neigh_modify    every 1 delay 0 check yes

thermo_style custom step temp press ke pe etotal
thermo_modify lost warn

thermo 1000

dump coords all atom 1 dump.xyz.dat

dump_modify coords scale no

restart 10000 restart.dat

timestep 1.00
run 20000

write_restart FINALCONFIG

"""