import numpy as np
import time
import pandas as pd
import random

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
        print('Incorrect format detected, rectifying the situation!')
        print('Assuming each coordinate is its own chain.')
        m = [(x,1) for x in range(len(no_dex))]
        multidex = pd.MultiIndex.from_tuples(m, names=['chain','atom'])
        no_dex = pd.DataFrame(data=no_dex.values, index=multidex)
    #maybe make this work for non homogenous systems later
    print('Assuming a homogenous system')
    no_mono = int(no_dex.index.get_level_values('atom').max()) 
    no_chains = int(no_dex.index.get_level_values('chain').max()) + 1
    print(f"chain number:\t{no_chains}\nmonomer number:\t{no_mono}")
    #rewrite this as a comprehension later
    tup_rep = []

    m1 = [[(i, 2*x) for x in range(no_mono)] for i in range(no_chains)]
    m1 = sum(m1, [])

    m2 = [[(i, 2*x+1) for x in range(no_mono)] for i in range(no_chains)]
    m2 = sum(m2, [])
    
    multidex_1 = pd.MultiIndex.from_tuples(m1, names=['chain', 'atom'])
    multidex_2 = pd.MultiIndex.from_tuples(m2, names=['chain', 'atom'])

    print(multidex_1)
    print(multidex_2)

    return multidex_1, multidex_2

def backmap(input_coordinates):

    multidex_1, multidex_2 = gen_multidex(input_coordinates)
    fudged = pd.DataFrame(fudge_position(input_coordinates))
    print(fudged)
    merged = pd.concat([
        pd.DataFrame(input_coordinates.values, multidex_1),
        pd.DataFrame(fudged.values, multidex_2)
    ])
    merged.sort_index(level='chain', inplace=True)
    merged.rename(columns={0:'x', 1:'y', 2:'z'}, inplace=True)
    print(merged)

    return merged


