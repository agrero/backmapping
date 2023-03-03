# module containing all functions relating to writing/reading files

import pandas as pd
import os

def write_lammps_input(filename, coordinates, bonds, angles, dihedrals, hi_lo, mass):
    """
    bad bad lazy programmer
    """
    parent = 'lammps_protocols'
    writename = os.path.join(parent, filename)
    
    no_atoms = len(coordinates)
    no_bonds = len(bonds)
    no_angles = len(angles)
    no_dihedrals = len(dihedrals)
    no_atom_types = 1
    no_bond_types = 1
    no_angle_types = 1 
    no_dihedral_types = 1

    xlo, xhi = hi_lo[0][0], hi_lo[0][1]
    ylo, yhi = hi_lo[1][0], hi_lo[1][1]
    zlo, zhi = hi_lo[2][0], hi_lo[2][1]
    header = f"""LAMMPS UA chain data file\n
{no_atoms} atoms
{no_bonds} bonds
{no_angles} angles
{no_dihedrals} dihedrals\n
{no_atom_types} atom types
{no_bond_types} bond types
{no_angle_types} angle types
{no_dihedral_types} dihedral types\n
{xlo} {xhi} xlo xhi
{ylo} {yhi} ylo yhi
{zlo} {zhi} zlo zhi\n
{'Masses'}\n
{'1'} {mass}
"""

    # theres a better way to do this
    with open(writename, 'w') as f:
        f.write(f'{header}\nAtoms\n\n')

    coordinates.round(4).to_csv(writename, sep=' ', mode='a', header=False)

    # checking to ensure these things actually exist
    # its science so its all fake but you get what i mean
    if no_bonds != 0:
        with open(writename, 'a') as f:
            f.write('\nBonds\n\n')
    
        bonds.to_csv(writename, sep=' ', mode='a', header=False)

    if no_angles != 0:
        with open(writename, 'a') as f:
            f.write('\nAngles\n\n')
    
        angles.to_csv(writename, sep=' ', mode='a', header=False)

    if no_dihedrals != 0:
        with open(writename, 'a') as f:
            f.write('\nDihedrals\n\n')
    
        dihedrals.to_csv(writename, sep=' ', mode='a', header=False)

def read_lammps(filename, components):
    """
    The first line of these out files aren't consistent so should just make the header the 
    initial row to the first section (Atoms # molecular) 

    Components are the components of the specific lammps file. IE Atoms, Bonds, Dihedrals.
    """
    parent = 'lammps_protocols'
    writename = os.path.join(parent, filename)
    with open(writename, 'r') as f:
        text = f.read()

        sections = text.split('\n')
        section_ndxs = [sections.index(i) for i in components]

        comp = []
        #header = sections[0:section_ndxs[0]]

        for ndx, i in enumerate(section_ndxs):
            try:
                slice = sections[i:section_ndxs[ndx+1]]
                comp.append(slice)
            except:
                slice = sections[i:]
                comp.append(slice)
        
        clean_sections = []
        for i in comp:
            ndxs = [1,-1]
            for j in ndxs:
                del i[j]
            clean_sections.append(i)
        #can change this later but will only return atoms for backmapping
        atoms = clean_sections[1][1:]
        cleaned_atoms = []
        #i should change it so i'm not writiing files with tabs but instead spaces
        for i in atoms:
            cleaned_atoms.append(i.split(' '))

        columns = ['atom', 'chain', 'atom type', 'x', 'y', 'z', 'xv', 'yv', 'zv']
        try:
            atom_frame = pd.DataFrame(data=cleaned_atoms, columns=columns)
        except:
            for i in atoms:
                cleaned_atoms.append(i.split('\t'))
            atom_frame = pd.DataFrame(data=cleaned_atoms, columns=columns)

        multidex = pd.MultiIndex.from_frame(atom_frame.iloc[:,:2])
        atom_frame.drop(['chain', 'atom','atom type','xv','yv','zv'], axis=1, inplace=True)
        print(atom_frame)

        atom_frame = pd.DataFrame(atom_frame.values, index=multidex, columns=['x','y','z']).astype(float)
        atom_frame.sort_index(level='chain', inplace=True)

        return atom_frame

def write_lammps_config(config_dict:dict, filename:str, restart=False, write_data=False, out_data='outdata', 
                        timesteps:list=[1], run_nos:list=[1]):
    """
    a function to write lammps configuration files using a dictionary as the main information source. 
    the individual inputs will be read in through the dictionary while the runnnig conditions will be generated
    through a function.

    config_dict: A dictionary consisting of LAMMPS input commands with the keys being the commands
        and values being the conditions.
    filename: a string filename to write the LAMMPS file to.
    restart: if true will read_restart, otherwise will read_data (i'm unsure if this really matters)
    write_data: this will write data post each run in order to make a little video!
    out_data: root naming pattern for write_data.
    timesteps: a list of timesteps for the lammps to run on,if not the same length as run_nos
        will extend itself by duplicating the last value until it is (unless if it is longer).
    run_nos: a list of number of runs to apply for each timestep, if not the same length as timesteps
        will extend itself by duplicating the last value until it is.
    """
    parent = 'lammps_protocols'
    writename = os.path.join(parent, filename)
    # checking whether or not to read restart or data
    if restart:
        del config_dict['read_data']
    else:
        try:
            del config_dict['read_restart']
        except:
            print('something was supposed to be here but aint, moving on')
    # if write_data is true will add that to the workflow
    # might wanna put this lower/after the run component
    if write_data:
        config_dict['write_data'] = f'write_data {out_data}'

    with open(writename, 'w') as f:
        for i in config_dict:
            f.write(f'{i} {config_dict[i]}\n')
        for ndx, i  in enumerate(timesteps):
            if ndx > 0:
                minimize = config_dict['minimize']
                f.write(f'minimize {minimize}\n')
            
            f.write(f'timestep {i}\n')
            f.write(f'run {run_nos[ndx]}\n')
        f.write('write_restart FINALCONFIG')

def write_backmapping_protocol(filename='protocol.sh', head_dict:dict=None, 
                               modules:list=['slurm', 'racs-spack', 'python3'],
                               environment:str="spack load /sfoi75e",
                               threads:int=1, no_iter:int=1, lammps_call:str='run'):
    """
    fill me in later lazy man
    """
    parent = 'lammps_protocols'
    writename = os.path.join(parent, filename)
    with open(writename, 'w') as f:
        f.write('#!/bin/bash\n')

    if head_dict != None:
        with open(writename, 'a') as f:
            for i in head_dict:
                if i == 'mail-type':
                    f.write(f'#SBATCH --{i}={head_dict[i][0]}\n')
                    f.write(f'#SBATCH --{i}={head_dict[i][1]}\n')
                    continue
                f.write(f'#SBATCH --{i}={head_dict[i]}\n')

    if modules != None:
        with open(writename, 'a') as f:
            f.write('BD=$PWD\ncd $BD\n')
            for module in modules:
                f.write(f'module load {module}\n')
    with open(writename, 'a') as f:
        for i in range(no_iter):
            # change this later to be able to make the input and outputfiles whatever you would want
            # might also  wanna make this less repetitive
            f.write(f'lmp < in.PEsoft > OUTsoft\n')
            f.write('cp FINALCONFIG RESTART\n')
            f.write(f'lmp < in.PEequil > OUTequil\n')
            f.write('cp FINALCONFIG RESTART\n')
            f.write('lmp read_restart FINALCONFIG\nlmp write_data outdata\n')
            f.write('python3 backmap.py\n')
            f.write(f'cp lammps-AT-config backmap-round-{i}\n')
            f.write('exit')