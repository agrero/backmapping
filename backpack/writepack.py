# module containing all functions relating to writing/reading files

import pandas as pd
import os
import mapback as mb

def write_lammps_input(filename, coordinates, hi_lo, mass):
    """
    bad bad lazy programmer
    """
    parent = 'lammps_protocols'
    writename = os.path.join(parent, filename)

    bonds = mb.make_bonds(coordinates)
    angles = mb.make_angles(coordinates)
    dihedrals = mb.make_dihedrals(coordinates)
    
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

def read_lammps_2(filepath):
    """
    We need to test to make sure that 'Atoms' is the beginning section before each 
    set of coordinates, honestly reading in the header isn't that important
    since we only need the coordinates and Masses

    ok cool so annotate this and do the mass thing
    """

    with open(filepath, 'r') as f:
        text = f.read()
        sections = text.split('\n')

        # pull the header from the text
        header_end_ndx = 0
        for ndx, i in enumerate(sections):
            if 'Atoms' in i:
                header_end_ndx = ndx
                break
        header = sections[0:header_end_ndx]

        # take the header and get the indices to slice the specific sections
        header_section_ndxs = []
        for ndx, i in enumerate(header):
            if len(header_section_ndxs) == 0:
                if i == '':
                    header_section_ndxs.append((0, ndx))
                    continue
            if i == '':
                header_section_ndxs.append((header_section_ndxs[-1][-1],ndx))

        header_sections = []
        for ndx, i in enumerate(header_section_ndxs):
            if ndx == 0:
                header_sections.append((header[i[0]:i[1]]))
            else:   
                header_sections.append((header[i[0]+1:i[1]]))

        summation = sum(header_sections[-2:], [])
        fixed_header_section = header_sections[:-2]
        fixed_header_section.append(summation)

        atom_no = int(fixed_header_section[1][0].split(' ')[0])

        atoms = sections[header_end_ndx+2:header_end_ndx+atom_no+2]
        atom_split = [i.split(' ') for i in atoms]

        columns = ['atom', 'chain', 'atom type', 'x', 'y', 'z', 'xv', 'yv', 'zv']
        if len(atom_split[0]) == 6:
            columns = columns[0:6]
            drop_columns = ['atom type']
            atom_frame = pd.DataFrame(data=atom_split, columns=columns)
            atom_frame = atom_frame.set_index(['chain', 'atom'])
            atom_frame.drop(drop_columns, axis=1, inplace=True)
        else:    
            drop_columns = ['atom type', 'xv', 'yv', 'zv']
            atom_frame = pd.DataFrame(data=atom_split, columns=columns)
            atom_frame = atom_frame.set_index(['chain', 'atom'])
            atom_frame.drop(drop_columns, axis=1, inplace=True)

        return atom_frame.astype(float)

#change me to write what's in line with what is working on Talapas
def write_backmapping_protocol(filename='protocol.sh',
                               head_dict:dict=None, 
                               threads=5,
                               modules:list=['slurm', 'racs-spack', 'python3'],
                               environment:str="spack load /sfoi75e",
                               no_iter:int=1):
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

def read_lammps_header(filepath):
    """
    reads a lammps header in order to ascertain which components are in the simulation, IE atoms,
    bonds, etc.

    filepath: the filepath to the lammps configuration/coordinates
    """
    components = []
    with open(filepath, 'r') as f:
        #may want to speed test which is faster, the short read or the just reading the first 150 bits
        text = f.read(150)
        split = text.split('\n')
        #this will give us the components of our lammps file
        split_dex = []
        for ndx, i in enumerate(split):
            if len(split_dex) == 2:
                break
            elif i == '':
                split_dex.append(ndx)


        components.append(split[(split_dex[0]+1):split_dex[1]])
    cleaned_components = []
    for i in components[0]:
        split = i.split(' ')
        cleaned_components.append(split[1])
            
    return cleaned_components


#should add a section in here to output movies! in line with what we have on talapas
def write_lammps_config(config_dict:dict, 
                        filename:str, 
                        initial_input:str='lammps-AT-config', 
                        restart=False, 
                        write_data=False, 
                        out_data='outdata', 
                        timesteps:list=[1], 
                        run_nos:list=[1]):
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
    # get the components in the lammps file
    parent = 'lammps_protocols'
    try:
        components = read_lammps_header(initial_input)
    except:
        components = read_lammps_header(os.path.join(parent, initial_input))
    #check to see out of all the lammps components which are not currently present in our input
    not_present = []
    for i in ['atoms', 'bonds', 'angles', 'dihedrals']:
        if not i in components:
            not_present.append(i)
    # get a list of commands to be parsed through
    lammps_commands = [i for i in config_dict.keys()]
    # parse lammps commands dict to see if there are commands for items not present
    for i in lammps_commands:
        for j in not_present:
            if j[0:-1] in i:
                del config_dict[i]


    parent = 'lammps_protocols'
    writename = os.path.join(parent, filename)

    if restart:
        del config_dict['read_data']
    else:
        try:
            del config_dict['read_restart']
        except:
            print('something was might have supposed to be here but aint, moving on')
    # if write_data is true will add that to the workflow
    # might wanna put this lower/after the run component

    if write_data:
        config_dict['write_data'] = f'{out_data}'
    
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