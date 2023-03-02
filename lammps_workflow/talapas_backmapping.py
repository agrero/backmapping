import numpy as np
import pandas as pd
import random
import time
import os

def fudge_position(data, radius=1):
    """
    data: this is a dataframe of xyz coordinates
    radius: the distance range of which to 'fudge' the coordinates

    duplicates the input data frame and 'fudges' each position at a distance between
    the -radius and +radius. 

    Returns a 'fudged' dataframe of the original input dataframe.
    """
    print('fudging location values')
    t0 = time.process_time_ns()
    fudgers = [[random.uniform(-radius, radius) for i in range(3)] for i in range(len(data))]

    fudgers = pd.DataFrame(fudgers, columns=['x','y','z'])
    data_xyz = data[['x','y','z']]

    new_data = data_xyz.values + fudgers
    t1 = time.process_time_ns()
    print(f'fudged locations in: {t1-t0}ns\n')

    return new_data

def gen_multidex(no_dex):
    """
    no_dex: a dataframe minimally containing x, y, and z coordinates

    if the dataframe is not formatted with a chain and atom multidex it will
    automatically generate said starting multidex. 

    backmaps each atom (of a homogenous system) by effectively doubling the length
    of the given dataframe. it will then generate two lists of indices
    of atoms n and n+1. it will maintain the chain # for each of the atoms.

    returns two multidexes as described earlier.
    """

    #correctly formatting indices
    print('indexing atoms')
    t0 = time.process_time_ns()
    if no_dex.index.names != ['atom', 'chain']:
        print('incorrect format detected, rectifying the situation!')
        print('assuming each coordinate is its own chain.')
        m = [(x,1) for x in range(len(no_dex))]
        multidex = pd.MultiIndex.from_tuples(m, names=['chain','atom'])
        no_dex = pd.DataFrame(data=no_dex.values, index=multidex)
        
    #maybe make this work for non homogenous systems later
    print('assuming a homogenous system')
    no_mono = len(no_dex.index.get_level_values('atom').unique())
    no_chains = len(no_dex.index.get_level_values('chain').unique())
    chain_length = int(no_mono / no_chains)
    print(f"chain number:\t{no_chains}\nmonomer number:\t{no_mono}")

    chain = sum([[x+1] * chain_length for x in range(0, no_chains)], [])
    atom_1 = [1 + 2*x for x in range(no_mono)]
    atom_2 = [1 + 2*x+1 for x in range(no_mono)]

    m1 = [(chain[i], atom_1[i]) for i in range(0, len(atom_1))]
    m2 = [(chain[i], atom_2[i]) for i in range(0, len(atom_2))]

    multidex_1 = pd.MultiIndex.from_tuples(m1, names=['chain', 'atom'])
    multidex_2 = pd.MultiIndex.from_tuples(m2, names=['chain', 'atom'])
    t1 = time.process_time_ns()
    print(f'multidex assembled in: {t1-t0}ns\n')
    
    return multidex_1, multidex_2   

def make_bonds(coordinates):
    """
    coordinates: a dataframe containing x, y, and z coordinates. indexed with a 
    chain and atom multidex.

    outputs a dataframe formatted as: index | type | ai | aj
    with ai and aj being the initial and final atoms in the bond
    """
    #index | type | ai | aj
    #retrieve the number of chains and how long each chain is
    print('assembling bonds')
    t0=time.process_time_ns()
    no_mono = len(coordinates.index.get_level_values('atom').unique())
    no_chains = len(coordinates.index.get_level_values('chain').unique())
    chain_length = int(no_mono/no_chains)

    ai = []
    for dn in range(no_chains):
        for n in range(chain_length-1):
            ai.append((1, n + dn*chain_length + 1, n + dn*chain_length + 2))

    bond_frame = pd.DataFrame(ai, columns=['type','ai','aj'], index=list(range(1, len(ai) + 1)))
    t1 = time.process_time_ns()
    print(f'bonds assembled in: {t1-t0}ns\n')
    return bond_frame

def reconfig_frame(coordinates):
    """
    coordinates: an x, y, and z coordinates dataframe

    reconfigures a frame so that the indexing is more in line with a lammps input 
    file by changing the index position and adding an atom type column.
    """
    print('reconfiguring frame')
    t0=time.process_time_ns()
    df = coordinates
    df = df.swaplevel(0, 1, axis=0)
    typ_array = [1] * (len(df))
    df.insert(0, 'type', typ_array)
    t1=time.process_time_ns()
    print(f'frame reconfigured in: {t1-t0}ns\n')
    return df

def make_angles(coordinates):
    """
    coordinates: a dataframe containing x, y, and z coordinates. indexed with a 
    chain and atom multidex.

    outputs a dataframe formatted as: index | type | ai | aj | ak
    with ai and aj being the initial and final atoms in the bond
        """
    #index | type | ai | aj
    #retrieve the number of chains and how long each chain is
    print('assembling angles')
    t0 = time.process_time_ns()
    no_mono = len(coordinates.index.get_level_values('atom').unique())
    no_chains = len(coordinates.index.get_level_values('chain').unique())
    chain_length = int(no_mono/no_chains)

    ai = []
    for dn in range(no_chains):
        for n in range(chain_length-2):
            ai.append((1, n + dn*chain_length + 1, n + dn*chain_length + 2,n + dn*chain_length + 3))
    t1 = time.process_time_ns()
    
    angle_frame = pd.DataFrame(ai, columns=['type','ai','aj','ak'], index=list(range(1, len(ai) + 1)))

    print(f'angles assembled in: {t1-t0}ns\n')

    return angle_frame



def make_dihedrals(coordinates):
    #index | type | ai | aj | ak | al
    """
    wow this looks familiar
    """
    print('assembling dihedrals')
    t0=time.process_time_ns()
    no_mono = len(coordinates.index.get_level_values('atom').unique())
    no_chains = len(coordinates.index.get_level_values('chain').unique())
    chain_length = int(no_mono/no_chains)

    ai = []
    for dn in range(no_chains):
        for n in range(chain_length-3):
            ai.append((1, n + dn*chain_length + 1, n + dn*chain_length + 2,n + dn*chain_length + 3, n + dn*chain_length + 4))
    dihedral_frame = pd.DataFrame(ai, columns=['type','ai','aj','ak', 'al'])
    t1=time.process_time_ns()
    print(f'dihedrals assembled in: {t1-t0}ns\n')

    return dihedral_frame
    


def backmap(input_coordinates):

    print('beginning backmapping\n')
    t0=time.process_time_ns()
    multidex_1, multidex_2 = gen_multidex(input_coordinates)
    
    fudged = pd.DataFrame(fudge_position(input_coordinates))
  

    merged = pd.concat([
        pd.DataFrame(input_coordinates.values, multidex_1),
        pd.DataFrame(fudged.values, multidex_2),
    ])
    merged.sort_index(level='chain', inplace=True)
    merged.rename(columns={0:'x', 1:'y', 2:'z'}, inplace=True)

    t1=time.process_time_ns()
    print(f'backmapping completed in: {t1-t0}ns\n')
    return merged

def write_lammps_input(filename, coordinates, bonds, angles, dihedrals, hi_lo, mass):
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

    #there may be a better way to do this but i just need it
    #to work
    with open(filename, 'w') as f:
        f.write(f'{header}\nAtoms\n\n')

    coordinates.round(4).to_csv(filename, sep=' ', mode='a', header=False)

    with open(filename, 'a') as f:
        f.write('\nBonds\n\n')
    
    bonds.to_csv(filename, sep=' ', mode='a', header=False)

    with open(filename, 'a') as f:
        f.write('\nAngles\n\n')
    
    angles.to_csv(filename, sep=' ', mode='a', header=False)

    with open(filename, 'a') as f:
        f.write('\nDihedrals\n\n')
    
    dihedrals.to_csv(filename, sep=' ', mode='a', header=False)

def read_lammps(filename, components):
    """
    The first line of these out files aren't consistent so should just make the header the 
    initial row to the first section (Atoms # molecular) 

    Components are the components of the specific lammps file. IE Atoms, Bonds, Dihedrals.
    """
    with open(filename, 'r') as f:
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
    