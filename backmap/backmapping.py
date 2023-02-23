import numpy as np
import pandas as pd
import random
import time

def fudge_position(data, radius=1):
    """
    the radius being divided by 100 will need to change
    """
    fudgers = [[random.uniform(-radius, radius) for i in range(3)] for i in range(len(data))]

    fudgers = pd.DataFrame(fudgers, columns=['x','y','z'])

    new_data = data.values + fudgers

    return new_data

def gen_multidex(no_dex):
    """input needs to be formatted with a chain and the atom count will be inferred
    should put something in there that
    """
    #i'm just going to duplicate the code below but should try to be more clever later
    #need to do this with a loop based on how many atoms are in a chain


    #correctly formatting indices
    if no_dex.index.names != ['chain', 'atom']:
        #print('Incorrect format detected, rectifying the situation!')
        #print('Assuming each coordinate is its own chain.')
        m = [(x,1) for x in range(len(no_dex))]
        multidex = pd.MultiIndex.from_tuples(m, names=['chain','atom'])
        no_dex = pd.DataFrame(data=no_dex.values, index=multidex)
    #maybe make this work for non homogenous systems later
    #print('Assuming a homogenous system')
    no_mono = len(no_dex.index.get_level_values('atom').unique())
    no_chains = len(no_dex.index.get_level_values('chain').unique())
    #print(f"chain number:\t{no_chains}\nmonomer number:\t{no_mono}")


    m1 = [[(i, 2*x) for x in range(no_mono)] for i in range(no_chains)]
    m1 = sum(m1, [])

    m2 = [[(i, 2*x+1) for x in range(no_mono)] for i in range(no_chains)]
    m2 = sum(m2, [])

    multidex_1 = pd.MultiIndex.from_tuples(m1, names=['chain', 'atom'])
    multidex_2 = pd.MultiIndex.from_tuples(m2, names=['chain', 'atom'])

    return multidex_1, multidex_2

def make_bonds(coordinates):
    #index | type | ai | aj
    #retrieve the number of chains and how long each chain is
    chains = coordinates.index.get_level_values('chain').drop_duplicates()
    at_p_chain = coordinates.index.get_level_values('atom').drop_duplicates().max()

    resized_chains = list(chains) * (at_p_chain)
    resized_chains.sort()
    ai = [i for i in range(0, at_p_chain)] * len(chains)
    aj = [i for i in range(1, at_p_chain + 1)] * len(chains)
    type_array = [1] * len(ai)

    bond_dict ={
        'type' : type_array,
        'ai': ai,
        'aj': aj
    }

    return pd.DataFrame(bond_dict)

def reconfig_frame(coordinates):
    df = coordinates.reset_index()
    df.drop(labels='atom', axis=1, inplace=True)
    typ_array = [1] * (len(df))
    df.insert(1, 'type', typ_array)

    return df

def make_angles(coordinates):
    #index | type | ai | aj | ak
    """
    stop being a lazy programmer

    on hold for the moment
    """
    at_p_chain = coordinates.index.get_level_values('atom').drop_duplicates().max()
    angle_atoms = [[i, i+1, i+2] for i in range(at_p_chain - 1)]

def make_dihedrals(coordinates):
    #index | type | ai | aj | ak | al
    pass
    


def backmap(input_coordinates):

    multidex_1, multidex_2 = gen_multidex(input_coordinates)
    
    fudged = pd.DataFrame(fudge_position(input_coordinates))
    typ_array = pd.Series(([1] * (len(fudged) * 2)))

    merged = pd.concat([
        pd.DataFrame(input_coordinates.values, multidex_1),
        pd.DataFrame(fudged.values, multidex_2),
    ])
    merged.sort_index(level='chain', inplace=True)
    merged.rename(columns={0:'x', 1:'y', 2:'z'}, inplace=True)

    return merged

def write_lammps_input(coordinates, hi_lo, mass):
    no_atoms = len(coordinates)
    no_atom_types = 1

    xlo, xhi = hi_lo[0][0], hi_lo[0][1]
    ylo, yhi = hi_lo[1][0], hi_lo[1][1]
    zlo, zhi = hi_lo[2][0], hi_lo[2][1]
    header = f"""
LAMMPS UA chain data file\n
{no_atoms} atoms\n
{no_atom_types} atom types\n
{xlo} {xhi} xlo xhi
{ylo} {yhi} ylo yhi
{zlo} {zhi} zlo zhi\n
{'Masses'}\n
{'1'} {mass}\n
"""

    #there may be a better way to do this but i just need it
    #to work
    with open('header_test.txt', 'w') as f:
        f.write(f'{header}\nAtoms\n')
        f.write('\n')

    coordinates.round(4).to_csv('header_test.txt', sep='\t', mode='a', header=False)
        



