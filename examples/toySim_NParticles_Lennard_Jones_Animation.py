import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from matplotlib.animation import FuncAnimation
import pandas as pd

plot_limit = 4

df = pd.read_excel('Data.xlsx')
print(df)

columns = list(df.columns.values)

#Find Number of Coordinates
coord_list = []
n=0
while n <= len(columns)-1:
    print(columns[n][0])
    if columns[n][0] == "x":
        coord_list.append(columns[n])
    n+=1
print(len(coord_list))

#Number of Atoms
N = len(coord_list)/3

fig, ax = plt.subplots(subplot_kw=dict(projection="3d"))

ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_zlim(-2, 2)

def update(theta):
    global point
    ax.clear()
    n=0
    time_coord_list_x = []
    time_coord_list_y = []
    time_coord_list_z = []
    while n <= N-1:
        time_coord_list_x.append(df.iloc[int(theta),2 + 3*n])
        time_coord_list_y.append(df.iloc[int(theta), 3 + 3 * n])
        time_coord_list_z.append(df.iloc[int(theta), 4 + 3 * n])
        ax.scatter(df.iloc[int(theta),2 + 3*n], df.iloc[int(theta),3 + 3*n], df.iloc[int(theta),4 + 3*n], color = 'black')
        # ax.quiver(df.iloc[int(theta), 2 + 3 * n], df.iloc[int(theta), 3 + 3 * n], df.iloc[int(theta), 4 + 3 * n], df.iloc[int(theta), 2*int(3*N) + 2 + 3 * n], df.iloc[int(theta), 2*int(3*N) + 3 + 3 * n], df.iloc[int(theta), 2*int(3*N) + 4 + 3 * n], color='red')
        n+=1
    ax.plot(time_coord_list_x, time_coord_list_y, time_coord_list_z, color = 'green')
    # ax.scatter(*get_coordinates(theta))
    ax.set_xlim(-plot_limit, plot_limit)
    ax.set_ylim(-plot_limit, plot_limit)
    ax.set_zlim(-plot_limit, plot_limit)
    ax.set_title("{}-Particles Lennard-Jones Simulation (T = {})".format(int(N),round(df["t"][theta],1)))
    # point.remove()
    # point = ax.scatter(*get_coordinates(theta))

ani = FuncAnimation(fig, update, frames=len(df["t"])-1, interval=50)
plt.show()
