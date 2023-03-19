# analysis module for backpack
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from backpack import mapback as mb
from backpack import writepack as wp

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
    distance.rename(index='distance', inplace=True)

    return distance

def get_bond_angles(coordinates):
    """intakes a set of coordinates in lammps format, generates sets of coordinates to calculate
    bond angles with them, and returns an array of said bond angles.

    coordinates: a set of xyz coordinates in lammps format (as below)
    atom | chain | x | y | z
    """

    # make angles 
    angles = mb.make_angles(coordinates)
    #print(angles)

    # get the indices of each angle
    no_chain = coordinates.reset_index(level='chain')
    no_chain.drop(labels=['chain'], axis=1, inplace=True)
    #print(no_chain)
    # may need to add a query here for ridiculous bond lengths 
    #   if the angle stuff is not fully apparent
    # remember to change this back to str when it works
    ai_ind = angles.loc[:,'ai'].values.astype(str)
    aj_ind = angles.loc[:,'aj'].values.astype(str)
    ak_ind = angles.loc[:,'ak'].values.astype(str)

    # for reference, ai is the set of xyz coordinates for atom i in the angle
    ai = no_chain[no_chain.index.isin(ai_ind)]
    ai.reset_index(drop=True, inplace=True)
    aj = no_chain[no_chain.index.isin(aj_ind)]
    aj.reset_index(drop=True, inplace=True)
    ak = no_chain[no_chain.index.isin(ak_ind)]
    ak.reset_index(drop=True, inplace=True)

    angles = []
    for ndx, i in enumerate(ai.values):
        #aj.loc[ndx].values
        ji = np.subtract(i,aj.loc[ndx].values)
        jk = np.subtract(ak.loc[ndx].values, aj.loc[ndx].values)

        cosine_angle = np.dot(ji, jk) / (np.linalg.norm(ji) * np.linalg.norm(jk))
        angle = np.arccos(np.around(cosine_angle,4))
        angles.append(np.degrees(angle))
    # maybe when there's less pressure redo it the way you did it before 
    # see how it compares, I think you'll probably just need to average the rows
    return pd.Series(angles)

def get_dihedral_angles(coordinates):
    """

    """
    dihedrals = mb.make_dihedrals(coordinates)

    no_chain = coordinates.reset_index(level='chain')
    no_chain.drop(labels=['chain'], axis=1, inplace=True)

    ai_ind = dihedrals.loc[:,'ai'].values.astype(str)
    aj_ind = dihedrals.loc[:,'aj'].values.astype(str)
    ak_ind = dihedrals.loc[:,'ak'].values.astype(str)
    al_ind = dihedrals.loc[:,'al'].values.astype(str)

    ai = no_chain[no_chain.index.isin(ai_ind)]
    ai.reset_index(drop=True, inplace=True)
    aj = no_chain[no_chain.index.isin(aj_ind)]
    aj.reset_index(drop=True, inplace=True)
    ak = no_chain[no_chain.index.isin(ak_ind)]
    ak.reset_index(drop=True, inplace=True)
    al = no_chain[no_chain.index.isin(al_ind)]
    al.reset_index(drop=True, inplace=True)

    dihedral_angles = []
    for ndx, i in enumerate(ai.values):
        #make this more in line with the other function
        p0 = i
        p1 = aj.loc[ndx].values
        p2 = ak.loc[ndx].values
        p3 = al.loc[ndx].values

        b0 = -1.0*(p1-p0)
        b1 = p2 - p1
        b2 = p3 - p2

        b0xb1 = np.cross(b0, b1)
        b1xb2 = np.cross(b2, b1)

        b0xb1_x_b1xb2 = np.cross(b0xb1, b1xb2)

        y = np.dot(b0xb1_x_b1xb2, b1)*(1.0/np.linalg.norm(b1))
        x = np.dot(b0xb1, b1xb2)

        dihedral_angles.append(np.degrees(np.arctan2(y,x)))

    return pd.DataFrame(dihedral_angles)
    
"""file = 'lammps-AT-config'
parent = 'lammps_protocols'
path = os.path.join(parent, file)
#data = wp.read_lammps_2(path)
initial_data = {
    'x' : [24.969, 24.044, 22.785, 21.951, 23.672, 22.881, 23.691, 22.557],
    'y' : [13.428, 12.661, 13.482, 13.670, 11.328, 10.326, 9.935, 9.096],
    'z' : [30.692, 29.808, 29.543, 30.431, 30.466, 29.620, 28.389, 30.459]
}
multi_dex = [(1,1),(2,1),(3,1),(4,1),(5,2),(6,2),(7,2),(8,2)]
mutlidex = pd.MultiIndex.from_tuples(multi_dex, names=['atom', 'chain'])
data = pd.DataFrame(initial_data, index=mutlidex)"""


def plot_bondlength(coordinates, filter_no=2, bin_width=.25):

    # get bond distances
    bond_distances = get_bond_len(coordinates)
    filtered = bond_distances.loc[lambda x : x < filter_no]

    # get dist
    bin_form = np.arange(0, max(filtered) + bin_width, bin_width)
    plt.hist(filtered, density=True, bins=bin_form)
    plt.xticks(bin_form)
    plt.ylabel('Probability')
    plt.xlabel('Bond Length')

    plt.show()

def plot_bondangle(coordinates, bin_width=1):

    # get bond angles
    bond_angle = get_bond_angles(coordinates)

    bin_form = np.arange(0, max(bond_angle) + bin_width, bin_width)

    plt.hist(bond_angle, density=True, bins=bin_form)
    plt.xticks(np.arange(0, max(bond_angle) + bin_width, bin_width * 2.))
    plt.ylabel('Probability')
    plt.xlabel('Bond Angle')

    plt.show()

def plot_dihedral(coordinates, bin_width=1):

    # get dihedral angles
    dihedral_angle = get_dihedral_angles(coordinates)

    bin_form = np.arange(-180, 180 + bin_width, bin_width)

    plt.hist(dihedral_angle, density=True, bins=bin_form)

    plt.xticks(np.arange(-180, 180 + bin_width, bin_width * 6))
    plt.ylabel('Probability')
    plt.xlabel('Dihedral Angle')

    plt.show()