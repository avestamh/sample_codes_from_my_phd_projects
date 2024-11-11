import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde

##_______read the file______________

data_file = np.genfromtxt('rmsd-rg-all.dat')

## __________generate your data_________

x = data_file[:,6] # rmsd_data

y = data_file[:,8] ## rg data

##___________ Calculate the point density__________________

xy = np.vstack([x,y])

z = gaussian_kde(xy)(xy)

##_____# Sort the points by density, so that the densest points are plotted last____

idx = z.argsort()

x, y, z = x[idx], y[idx], z[idx]

fig = plt.figure(figsize=(6,6))

grid = plt.GridSpec(4, 4, hspace=0.3, wspace=0.6)

main_ax = fig.add_subplot(grid[:-1, 1:])

y_hist = fig.add_subplot(grid[:-1, 0], xticklabels=[], sharey=main_ax)

x_hist = fig.add_subplot(grid[-1, 1:], yticklabels=[], sharex=main_ax)

##______plotting the density map_____

dens_map = main_ax.scatter(x, y, c=z, s=30, cmap='jet')

### histogram on the echa sides of the density map

x_hist.hist(x, 100, histtype='stepfilled',orientation='vertical', color='b')

x_hist.invert_yaxis()

y_hist.hist(y, 100, histtype='stepfilled', orientation='horizontal', color='b')

y_hist.invert_xaxis()

### set the labels

y_hist.set_ylabel(r'Rg ($\AA$)', fontsize=20)

x_hist.set_xlabel(r'RMSD ($\AA$)', fontsize=20)

##_________ adjust the colorbar posotions_________

# fig.subplots_adjust(right=0.7)

cbar_ax = fig.add_axes([0.92, 0.33, 0.03, 0.55])

plt.colorbar(dens_map, cax=cbar_ax)

img_name= 'rmsd-rg-dense-hist.png'

plt.savefig(img_name, dpi=300, bbox_inches='tight')

# plt.show()
