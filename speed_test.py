import numpy as np
import pandas as pd

import time
import os

import matplotlib.pyplot as plt

from backmap import backmap as bm

def generate_test_data(input_coordinates, no_iterations):
    os.chdir("C:\Coding Projects\Guenza Lab\\backmapping\csvs")
    backmapped = input_coordinates
    for i, val in enumerate(list(range(no_iterations))):

        backmapped = bm.backmap(backmapped)
        no_mono = len(backmapped.index.get_level_values('atom').unique())
        no_chains = len(backmapped.index.get_level_values('chain').unique())
        filename = f'backmap_{i}_mononum_{no_mono}_chainnum_{no_chains}.csv'
        backmapped.to_csv(filename)

def speed_test(directory_path, no_iterations):
    dir_list = os.listdir(directory_path)
    os.chdir(directory_path)
    x = []
    y = []
    
    for i in dir_list:
        to_backmap = pd.read_csv(i, index_col=['chain', 'atom'])
        no_mono = len(to_backmap.index.get_level_values('atom').unique())
        no_chains = len(to_backmap.index.get_level_values('chain').unique())
        x.append(no_mono * no_chains)

        timings = []
        for j in range(no_iterations):
            t0 = time.time()
            bm.backmap(to_backmap)
            t1 = time.time()
            tf = t1-t0
            timings.append(tf)
        ave = np.average(np.array(timings))
        y.append(ave)
        timings.clear()

    return x , y


#data = {
#    'x':[1.0,1.05,1.2],
#    'y':[1.0,1.05,1.2],
#    'z':[1.0,1.05,1.2]
#}

#init_coordinate = pd.DataFrame(data, columns=data.keys())
#generate_test_data(init_coordinate,50)

dir_path = "C:\Coding Projects\Guenza Lab\\backmapping\csvs"
x,y = speed_test(dir_path, 1)

plt.style.use('_mpl-gallery')

fig, ax = plt.subplots(figsize = (9,9))
ax.margins(x=0.5, y=0.5)
ax.scatter(x, y, s=60, edgecolors='k')

plt.xlabel('Total Number of Monomers')
plt.ylabel('Time (s)')
plt.show()



"""
for i in range(500):
    t0 = time.time()
    backmap(init_coordinate)
    t1 = time.time()
    tf = t1-t0
    timings.append(tf)

data = np.array(timings)
ave = np.average(data)
print(ave)"""