# analysis module for backpack
import os

import mapback as mb
import writepack as wp
import numpy as np

def get_bond_len(coordinates):
    """intakes a set of coordinates in lammps format, creates bonds for all of them,
    and returns an array of said bond distances.
    
    coordinates: a set of xyz coordinates in lammps format 
    atom | chain | x | y | z
    """
    # make bonds 
    bonds = mb.make_bonds(coordinates)

    # get the indices of each bond 
    no_chain = coordinates.reset_index(level='chain')
    no_chain.drop(labels=['chain'], axis=1, inplace=True)
    ai = bonds.loc[:,'ai'].values.astype(str)
    aj = bonds.loc[:,'aj'].values.astype(str)

    # from previously retrieved indices, get the atom coordinates 
    xyz_1 = no_chain[no_chain.index.isin(ai)]
    xyz_1.reset_index(drop=True, inplace=True)
    xyz_2 = no_chain[no_chain.index.isin(aj)]
    xyz_2.reset_index(drop=True, inplace=True)

    # can rewrite as a loop over x, y, z
    #
    power = np.repeat([2], len(xyz_1))

    x_diff = np.power(np.abs(np.subtract(xyz_1['x'], xyz_2['x'])), power)
    y_diff = np.power(np.abs(np.subtract(xyz_1['y'], xyz_2['y'])), power)
    z_diff = np.power(np.abs(np.subtract(xyz_1['z'], xyz_2['z'])), power)

    distance = np.sqrt(np.abs(np.subtract(np.subtract(x_diff, y_diff), z_diff)))
    return distance

def get_bond_angles(coordinates):
    """intakes a set of coordinates in lammps format, generates sets of coordinates to calculate
    bond angles with them, and returns an array of said bond angles.

    coordinates: a set of xyz coordinates in lammps format (as below)
    atom | chain | x | y | z
    """

    # make angles 
    angles = mb.make_angles(coordinates)


    # get the indices of each angle
    no_chain = coordinates.reset_index(level='chain')
    no_chain.drop(labels=['chain'], axis=1, inplace=True)
    # may need to add a query here for ridiculous bond lengths 
    #   if the angle stuff is not fully apparent
    ai_ind = angles.loc[:,'ai'].values.astype(str)
    aj_ind = angles.loc[:,'aj'].values.astype(str)
    ak_ind = angles.loc[:,'ak'].values.astype(str)

    ai = no_chain[no_chain.index.isin(ai_ind)]
    ai.reset_index(drop=True, inplace=True)
    aj = no_chain[no_chain.index.isin(aj_ind)]
    aj.reset_index(drop=True, inplace=True)
    ak = no_chain[no_chain.index.isin(ak_ind)]
    ak.reset_index(drop=True, inplace=True)
    
    # do math
    ba = ai - aj
    bc = ak - aj

    cos_angle = (ba.T * bc.T) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cos_angle)

    return np.degrees(angle)

 
file = 'lammps-AT-config'
parent = 'lammps_protocols'
path = os.path.join(parent, file)
#data = wp.read_lammps_2(path)
data = mb.generate_test_data(3, 23.4711697)
print(data)
data = mb.backmap(data, 23.4711697)
bond_len = get_bond_len(data)
bond_angle = get_bond_angles(data)
print(bond_angle.head(10))
#print(np.average(bond_angle))

