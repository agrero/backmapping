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
    print('fudging location values')
    t0 = time.time()
    fudgers = [[random.uniform(-radius, radius) for i in range(3)] for i in range(len(data))]

    fudgers = pd.DataFrame(fudgers, columns=['x','y','z'])
    data_xyz = data[['x','y','z']]

    new_data = data_xyz.values + fudgers
    t1 = time.time()
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

    # correctly formatting indices
    print('indexing atoms')
    t0 = time.time()
    print(no_dex.index.names)
    # this shit is kinda broken
    if no_dex.index.names[1] != 'atom' or no_dex.index.names[0] != 'chain':
        print('incorrect format detected. processing. PROCESSING!')
        print('assuming each coordinate is its own chain.')
        m = [(x,x) for x in range(len(no_dex))]
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
    t1 = time.time()
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
    t0=time.time()
    no_mono = len(coordinates.index.get_level_values('atom').unique())
    no_chains = len(coordinates.index.get_level_values('chain').unique())
    chain_length = int(no_mono/no_chains)

    ai = []
    for dn in range(no_chains):
        for n in range(chain_length-1):
            ai.append((1, n + dn*chain_length + 1, n + dn*chain_length + 2))

    bond_frame = pd.DataFrame(ai, columns=['type','ai','aj'], index=list(range(1, len(ai) + 1)))
    t1 = time.time()
    print(f'bonds assembled in: {t1-t0}ns\n')
    return bond_frame

def reconfig_frame(coordinates):
    """
    coordinates: an x, y, and z coordinates dataframe

    reconfigures a frame so that the indexing is more in line with a lammps input 
    file by changing the index position and adding an atom type column.
    """
    print('reconfiguring frame')
    t0=time.time()
    df = coordinates
    df = df.swaplevel(0, 1, axis=0)
    typ_array = [1] * (len(df))
    df.insert(0, 'type', typ_array)
    t1=time.time()
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
    t0 = time.time()
    no_mono = len(coordinates.index.get_level_values('atom').unique())
    no_chains = len(coordinates.index.get_level_values('chain').unique())
    chain_length = int(no_mono/no_chains)

    ai = []
    for dn in range(no_chains):
        for n in range(chain_length-2):
            ai.append((1, n + dn*chain_length + 1, n + dn*chain_length + 2,n + dn*chain_length + 3))
    t1 = time.time()
    
    angle_frame = pd.DataFrame(ai, columns=['type','ai','aj','ak'], index=list(range(1, len(ai) + 1)))

    print(f'angles assembled in: {t1-t0}ns\n')

    return angle_frame

def make_dihedrals(coordinates):
    #index | type | ai | aj | ak | al
    """
    wow this looks familiar
    """
    print('assembling dihedrals')
    t0=time.time()
    no_mono = len(coordinates.index.get_level_values('atom').unique())
    no_chains = len(coordinates.index.get_level_values('chain').unique())
    chain_length = int(no_mono/no_chains)

    ai = []
    for dn in range(no_chains):
        for n in range(chain_length-3):
            ai.append((1, n + dn*chain_length + 1, n + dn*chain_length + 2,n + dn*chain_length + 3, n + dn*chain_length + 4))
    dihedral_frame = pd.DataFrame(ai, columns=['type','ai','aj','ak', 'al'])
    t1=time.time()
    print(f'dihedrals assembled in: {t1-t0}ns\n')

    return dihedral_frame


def check_validity(dataframe, box_bounds):
    query = ''
    xyz= ['x', 'y', 'z']
    for i in xyz:
        if i == 'z':
            query += f'{i} > {box_bounds} or {i} < -{box_bounds}'
        else:    
            query += f'{i} > {box_bounds} or {i} < -{box_bounds} or '
    # this might be the fastest way to do it, it'll get fixed when the system equils 
    print(len(dataframe.query(query)))
    for i in xyz:
        dataframe.loc[dataframe[f'{i}'] > box_bounds, f'{i}'] = box_bounds
        dataframe.loc[dataframe[f'{i}'] < -box_bounds, f'{i}'] = -box_bounds

    return dataframe
    
def backmap(input_coordinates, box_bounds):

    print('beginning backmapping\n')
    t0=time.time()
    multidex_1, multidex_2 = gen_multidex(input_coordinates)

    fudged = pd.DataFrame(fudge_position(input_coordinates))

    merged = pd.concat([
        pd.DataFrame(input_coordinates.values, multidex_1),
        pd.DataFrame(fudged.values, multidex_2),
    ])
    #problem occurs before here
    merged.sort_index(level='chain', inplace=True)
    merged.rename(columns={0:'x', 1:'y', 2:'z'}, inplace=True)
    print('checking boundary conditions')
    merged = check_validity(merged, box_bounds)

    t1=time.time()
    print(f'backmapping completed in: {t1-t0}ns\n')
    return merged

def generate_test_data(no_monomers=10, bounds=10):
    xyz = [[random.uniform(-bounds, bounds) for i in range(3)] for i in range(no_monomers)]
    coordinates = pd.DataFrame(xyz, columns=['x','y','z'])
    return backmap(coordinates)
