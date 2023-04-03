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

    if type(list(no_chain.index.values)[0]) == np.dtype('int64'):
        ai = bonds.loc[:,'ai'].values.astype(int)
        aj = bonds.loc[:,'aj'].values.astype(int)
    elif type(list(no_chain.index.values)[0]) == str:
        ai = bonds.loc[:,'ai'].values.astype(str)
        aj = bonds.loc[:,'aj'].values.astype(str)
    else:
        raise ValueError('expected np int64 or string for indices')
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

    # may need to add a query here for ridiculous bond lengths 
    #   if the angle stuff is not fully apparent
    # remember to change this back to str when it works
    if type(list(no_chain.index.values)[0]) == np.dtype('int64'):
        ai_ind = angles.loc[:,'ai'].values.astype(int)
        aj_ind = angles.loc[:,'aj'].values.astype(int)
        ak_ind = angles.loc[:,'ak'].values.astype(int)
    elif type(list(no_chain.index.values)[0]) == str:        
        ai_ind = angles.loc[:,'ai'].values.astype(str)
        aj_ind = angles.loc[:,'aj'].values.astype(str)
        ak_ind = angles.loc[:,'ak'].values.astype(str)
    else:
        raise ValueError('expected np int64 or string for indices')      

    # for reference, ai is the set of xyz coordinates for atom i in the angle
    ai = no_chain[no_chain.index.isin(ai_ind)]
    ai.reset_index(drop=True, inplace=True)
    aj = no_chain[no_chain.index.isin(aj_ind)]
    aj.reset_index(drop=True, inplace=True)
    ak = no_chain[no_chain.index.isin(ak_ind)]
    ak.reset_index(drop=True, inplace=True)

    # calculating the bond angles
    angles = []
    for ndx, i in enumerate(ai.values):
        ji = np.subtract(i,aj.loc[ndx].values)
        jk = np.subtract(ak.loc[ndx].values, aj.loc[ndx].values)

        cosine_angle = np.dot(ji, jk) / (np.linalg.norm(ji) * np.linalg.norm(jk))
        angle = np.arccos(np.around(cosine_angle,4))
        angles.append(np.degrees(angle))
    # maybe when there's less pressure redo it the way you did it before 
    # see how it compares, I think you'll probably just need to average the rows
    return pd.Series(angles)

def get_dihedral_angles(coordinates):
    """intakes a set of coordinates in lammps format, generates sets of coordinates to calculate
    dihedral angles with them, and returns an array of said bond angles.

    coordinates: a set of xyz coordinates in lammps format (as below)
    atom | chain | x | y | z 
    """
    # get dihedral angles 
    dihedrals = mb.make_dihedrals(coordinates)

    no_chain = coordinates.reset_index(level='chain')
    no_chain.drop(labels=['chain'], axis=1, inplace=True)
    if type(list(no_chain.index.values)[0]) == np.dtype('int64'):
        ai_ind = dihedrals.loc[:,'ai'].values.astype(int)
        aj_ind = dihedrals.loc[:,'aj'].values.astype(int)
        ak_ind = dihedrals.loc[:,'ak'].values.astype(int)
        al_ind = dihedrals.loc[:,'al'].values.astype(int)
    elif type(list(no_chain.index.values)[0]) == str:
        ai_ind = dihedrals.loc[:,'ai'].values.astype(str)
        aj_ind = dihedrals.loc[:,'aj'].values.astype(str)
        ak_ind = dihedrals.loc[:,'ak'].values.astype(str)
        al_ind = dihedrals.loc[:,'al'].values.astype(str)
    else:
        raise ValueError('expected np int64 or string for indices')  

    ai = no_chain[no_chain.index.isin(ai_ind)]
    ai.reset_index(drop=True, inplace=True)
    aj = no_chain[no_chain.index.isin(aj_ind)]
    aj.reset_index(drop=True, inplace=True)
    ak = no_chain[no_chain.index.isin(ak_ind)]
    ak.reset_index(drop=True, inplace=True)
    al = no_chain[no_chain.index.isin(al_ind)]
    al.reset_index(drop=True, inplace=True)

    # calculating dihedral angles
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

def plot_bondlength(coordinates, filter_no=2, bin_width=.25):
    """
    intakes the output of the bondlength function and plots the distribution of bond lengths 
    as a count distribution.

    coordinates: dataframe of bond lengths as explained above
    filter_no: cutoff for bond length
    bin_width: range of each bin the histogram uses.
    """
    # get bond distances
    bond_distances = get_bond_len(coordinates)
    filtered = bond_distances.loc[lambda x : x < filter_no]

    # get dist
    bin_form = np.arange(0, max(filtered) + bin_width, bin_width)

    y, binedges = np.histogram(filtered, bins=bin_form)

    plt.hist(filtered, bins=bin_form, edgecolor='orange')

    bincenters = 0.5 * (binedges[1:] + binedges[:-1])
    plt.plot(bincenters, y, '-', c='black')
    plt.xticks(np.arange(0, max(filtered) + bin_width, bin_width * xtick_dist))
    plt.ylabel('Percentage')
    plt.xlabel('Bond Length')

    plt.show()

def plot_bondangle(coordinates, bin_width=1, xtick_dist=5):
    """
    inputs the output of the bondangle function and plots it as histogram of indivdual counts

    coordinates: bond angles as stated above
    bin_width: range of each bin the histogram uses.
    """
    # get bond angles
    bond_angle = get_bond_angles(coordinates)

    bin_form = np.arange(min(bond_angle), max(bond_angle) + bin_width, bin_width)

    plt.hist(bond_angle, density=True, bins=bin_form)
    plt.xticks(np.arange(min(bond_angle), max(bond_angle) + bin_width, bin_width * xtick_dist))
    plt.ylabel('Percentage')
    plt.xlabel('Bond Angle')

    plt.show()

def plot_dihedral(coordinates, bin_width=1):
    """
    inputs the output of the bondangle function and plots it as histogram of indivdual counts

    coordinates: dihedral angles as stated above
    bin_width: range of each bin the histogram uses.
    """
    # get dihedral angles
    dihedral_angle = get_dihedral_angles(coordinates)

    bin_form = np.arange(-180, 180 + bin_width, bin_width)

    plt.hist(dihedral_angle, density=True, bins=bin_form)

    plt.xticks(np.arange(-180, 180 + bin_width, bin_width * xtick_dist))
    plt.ylabel('Percentage')
    plt.xlabel('Dihedral Angle')

    plt.show()