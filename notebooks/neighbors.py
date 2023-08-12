import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

sdf = pd.read_csv('../Data/tree_algorithm/set_particles.csv')
part = pd.read_csv('../Data/tree_algorithm/set_neighbors.csv')
particle = part["p"][0]
neighbors = part["p"][1:]

x_ip = sdf["x"]-sdf.iloc[particle]["x"]
y_ip = sdf["y"]-sdf.iloc[particle]["y"]
z_ip = sdf["z"]-sdf.iloc[particle]["z"]
r_ip = np.sqrt(x_ip*x_ip + y_ip*y_ip + z_ip*z_ip)
radius = 2*sdf.iloc[particle]["h"]

fig, ax = plt.subplots(figsize=(5,5))
ax = fig.add_subplot(projection='3d')

real_neigh = []
for ii in range(len(r_ip)):
    if r_ip[ii] <= radius:
        real_neigh.append(ii)
        ax.scatter([sdf['x'][ii]], [sdf['y'][ii]], [sdf['z'][ii]], color='blue', marker=".", s=1)

for ii in neighbors:
    ax.scatter([sdf['x'][ii]], [sdf['y'][ii]], [sdf['z'][ii]], color='r', marker=".", s=1)
    
ax.scatter(sdf.iloc[particle]["x"], sdf.iloc[particle]["y"], sdf.iloc[particle]["z"], marker="o", color="g",s=5)

r = radius
u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
x = r*np.cos(u) * np.sin(v) + sdf.iloc[particle]["x"] 
y = r*np.sin(u) * np.sin(v) + sdf.iloc[particle]["y"]
z = r*np.cos(v) + sdf.iloc[particle]["z"]
ax.plot_surface(x, y, z, color="r", alpha=0.2)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

ax.set_xlim3d(0,1)
ax.set_ylim3d(0.1)
ax.set_zlim3d(0.1)
ax.set_box_aspect([1,1,1])

plt.show()