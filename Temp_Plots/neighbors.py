import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

sdf = pd.read_csv('../Data/tree_algorithm/set_particles.csv')
part = pd.read_csv('../Data/tree_algorithm/set_neighbors.csv')
particle = part["p"][0]
neighbors = part["p"][1:]

fig, ax = plt.subplots(figsize=(5,5))
ax = fig.add_subplot(projection='3d')
#plt.plot(sdf['x'][neightbors], sdf['y'][particle], 'o', color='crimson', ms=5)

for ii in neighbors:
    ax.scatter([sdf['x'][ii]], [sdf['y'][ii]], [sdf['z'][ii]], color='blue', marker=".")

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("y")

ax.set_xlim(0,1)
ax.set_ylim(0.1)
ax.set_zlim(0.1)

plt.show()