import matplotlib.pyplot as plt
import sarracen as src
import numpy as np
import os

file            = "toy_star"                            # Test's name to be animated 
data_address    = "../Data/results/" + file + "/"       # Full address of results files
images_address  = "../Data/results/" + file + "/Images/"# Full address of results files
time            = np.loadtxt(data_address+"Time.txt")   # Time data
T               = len(time)                             # Time steps

hfact   = 1.2
dm      = 1./src.read_csv("../Data/initial_distribution/"+file+".csv").shape[0]
fmratei = 2
fmrateo = 25

# Find upper limit of density 
rhomax = 0.0
for ii in range(T):
    sdf = src.read_csv(data_address + str(ii) + ".csv")
    sdf.params={'mass':dm, 'hfact':hfact}
    sdf.calc_density()
    rhotemp = sdf["rho"].max()
    if rhotemp > rhomax:
        rhomax = rhotemp

os.system(f"mkdir {data_address}Images")

for ii in range(T):
    sdf = src.read_csv(data_address + str(ii) + ".csv")
    sdf.params={'mass':dm, 'hfact':hfact}
    sdf.calc_density()
    fig, ax = plt.subplots(figsize=(7,5))
    ax = sdf.render('rho', ax=ax, rotation=[0,0,0], xlim=[-0.7, 0.7], ylim=[-0.7, 0.7], vmin=0.0, vmax=rhomax)
    ax.set_facecolor('k')
    ax.set_xlim(-0.7,0.7)
    ax.set_ylim(-0.7,0.7)
    ax.set_xticks(np.linspace(-0.6, 0.6, 7))
    ax.set_yticks(np.linspace(-0.6, 0.6, 7))
    ax.text(-0.6, 0.6, r"$t=$"+f"{time[ii]:.2f}", color="w")
    plt.savefig(images_address + str(ii) + ".png")
    plt.close()

os.system(f"ffmpeg -framerate {fmratei} -i '{images_address}%d.png' -start_number 0 -c:v libx264 -r {fmrateo} {data_address}evolution.mp4")
os.system(f"rm -r {images_address}")