import pandas as pd
import random
import time

# with the new way of doing this all through python we 
# should really get rid of this and go back to the package method 

def fudge_position(data, radius=1):
    """
    data: this is a dataframe of xyz coordinates
    radius: the distance range of which to 'fudge' the coordinates

    duplicates the input data frame and 'fudges' each position at a distance between
    the -radius and +radius. 

    Returns a 'fudged' dataframe of the original input dataframe.
    """

    fudgers = [[random.uniform(-radius, radius) for i in range(3)] for i in range(len(data))]

    fudgers = pd.DataFrame(fudgers, columns=['x','y','z'])
    data_xyz = data[['x','y','z']]

    new_data = data_xyz.values + fudgers

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

    # correctly formatting indices
    if no_dex.index.names[1] != 'atom' or no_dex.index.names[0] != 'chain':
        print('incorrect format detected. processing. PROCESSING!')
        print('assuming each coordinate is its own chain.')
        m = [(x,x) for x in range(len(no_dex))]
        multidex = pd.MultiIndex.from_tuples(m, names=['chain','atom'])
        no_dex = pd.DataFrame(data=no_dex.values, index=multidex)

    no_mono = len(no_dex.index.get_level_values('atom').unique())
    no_chains = len(no_dex.index.get_level_values('chain').unique())
    chain_length = int(no_mono / no_chains)

    chain = sum([[x+1] * chain_length for x in range(0, no_chains)], [])
    atom_1 = [1 + 2*x for x in range(no_mono)]
    atom_2 = [1 + 2*x+1 for x in range(no_mono)]


    m1 = [(chain[i], atom_1[i]) for i in range(0, len(atom_1))]
    m2 = [(chain[i], atom_2[i]) for i in range(0, len(atom_2))]

    multidex_1 = pd.MultiIndex.from_tuples(m1, names=['chain', 'atom'])
    multidex_2 = pd.MultiIndex.from_tuples(m2, names=['chain', 'atom'])

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

    no_mono = len(coordinates.index.get_level_values('atom').unique())
    no_chains = len(coordinates.index.get_level_values('chain').unique())
    chain_length = int(no_mono/no_chains)

    ai = []
    for dn in range(no_chains):
        for n in range(chain_length-1):
            ai.append((1, n + dn*chain_length + 1, n + dn*chain_length + 2))

    bond_frame = pd.DataFrame(ai, columns=['type','ai','aj'], index=list(range(1, len(ai) + 1)))

    return bond_frame

def reconfig_frame(coordinates):
    """
    coordinates: an x, y, and z coordinates dataframe

    reconfigures a frame so that the indexing is more in line with a lammps input 
    file by changing the index position and adding an atom type column.
    """

    df = coordinates
    df = df.swaplevel(0, 1, axis=0)
    typ_array = [1] * (len(df))
    df.insert(0, 'type', typ_array)
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

    no_mono = len(coordinates.index.get_level_values('atom').unique())
    no_chains = len(coordinates.index.get_level_values('chain').unique())
    chain_length = int(no_mono/no_chains)

    ai = []
    for dn in range(no_chains):
        for n in range(chain_length-2):
            ai.append((1, n + dn*chain_length + 1, n + dn*chain_length + 2,n + dn*chain_length + 3))

    angle_frame = pd.DataFrame(ai, columns=['type','ai','aj','ak'], index=list(range(1, len(ai) + 1)))

    return angle_frame

def make_dihedrals(coordinates):
    #index | type | ai | aj | ak | al
    """
    from a given set of 3d coordinates in lammps format, generate a set of dihedral angles.

    coordinates: 3d coordiantes in lammps format
    """

    no_mono = len(coordinates.index.get_level_values('atom').unique())
    no_chains = len(coordinates.index.get_level_values('chain').unique())
    chain_length = int(no_mono/no_chains)

    ai = []
    for dn in range(no_chains):
        for n in range(chain_length-3):
            ai.append((1, n + dn*chain_length + 1, n + dn*chain_length + 2,n + dn*chain_length + 3, n + dn*chain_length + 4))
    dihedral_frame = pd.DataFrame(ai, columns=['type','ai','aj','ak', 'al'])

    return dihedral_frame


def check_validity(dataframe, box_bounds):
    """
    checks to see that all of the given atoms are within the bounds of the simulation box,
    if they are not it will move the atom just into the frame of the box.

    dataframe: xyz coordinates in lammps format
    box_bounds: the square bounds of the simulation box (IE one side length)
    """
    query = ''
    xyz= ['x', 'y', 'z']
    for i in xyz:
        if i == 'z':
            query += f'{i} > {box_bounds} or {i} < -{box_bounds}'
        else:    
            query += f'{i} > {box_bounds} or {i} < -{box_bounds} or '
    # this might be the fastest way to do it, it'll get fixed when the system equils 

    for i in xyz:
        dataframe.loc[dataframe[f'{i}'] > box_bounds, f'{i}'] = box_bounds
        dataframe.loc[dataframe[f'{i}'] < -box_bounds, f'{i}'] = -box_bounds

    return dataframe
    
def backmap(input_coordinates, box_bounds):
    """
    the backmapping workflow

    input_coordinates: xyz coordinates in lammps format
    box_bounds: the lenght of one side of the square simulation box
    """
    multidex_1, multidex_2 = gen_multidex(input_coordinates)

    fudged = pd.DataFrame(fudge_position(input_coordinates))

    merged = pd.concat([
        pd.DataFrame(input_coordinates.values, multidex_1),
        pd.DataFrame(fudged.values, multidex_2),
    ])

    merged.sort_index(level='chain', inplace=True)
    merged.rename(columns={0:'x', 1:'y', 2:'z'}, inplace=True)

    merged = check_validity(merged, box_bounds)

    return merged

def generate_test_data(no_monomers=10, bounds=10):
    """generates a random set of monomers within a given square box bounds.

    no_monomers: number of course grained monomers
    bounds: the length of an individual side of the simulation box    
    """
    xyz = [[random.uniform(-bounds, bounds) for i in range(3)] for i in range(no_monomers)]
    coordinates = pd.DataFrame(xyz, columns=['x','y','z'])
    multi_dex = gen_multidex(coordinates)
    coordinates.reindex(multi_dex)
    return coordinates
