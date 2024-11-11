import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde 

# from mpl_toolkits.axes_grid1 import make_axes_locatable

# 

cmaps = {1:'Greens', 2:'Oranges', 3:'Purples',4:'Greys' }

hist_col = {1:'Green', 2:'Orange', 3:'Purple',4:'Grey' }

# print(pos)

fig = plt.figure(figsize=(10,10))

grid = plt.GridSpec(4, 4, hspace=0.3, wspace=0.45)

#  

for i in range(1,5):

    data = np.genfromtxt('data%d.dat'%i)

    y = data[0:2000, 8]

    x = data[0:2000, -1]

# 

    xy = np.vstack([x,y])

    z = gaussian_kde(xy)(xy)

    idx = z.argsort()

    x, y, z = x[idx], y[idx], z[idx]

    main_ax = plt.subplot(grid[:-1, 1:])

    y_hist = plt.subplot(grid[:-1, 0], xticklabels=[], sharey=main_ax)

    x_hist = plt.subplot(grid[-1, 1:], yticklabels=[], sharex=main_ax)

    dense_map = main_ax.scatter(x, y, c=z, s=20, cmap=cmaps[i], edgecolor='')

    x_hist.hist(x, 100, histtype='stepfilled', orientation='vertical',color=hist_col[i], alpha=.3, density='normed')

    # x_hist.invert_yaxis()

    y_hist.hist(y, 100, histtype='stepfilled', orientation='horizontal', color=hist_col[i], alpha=.3, density='normed')

    # y_hist.invert_xaxis()

    y_hist.set_ylabel(r'Rg ($\AA$)', fontsize=26)

    x_hist.set_xlabel(r'$\theta (^{\circ})$', fontsize=26)

    main_ax.set_title('Density Plot', fontsize=26)

    pos = main_ax.get_position().get_points().flatten()

    print('data%d: '%i, pos)

    if i == 1 or i ==2:

        cabar_position = [(pos[2]+0.005 +(0.055*(i%2))), (pos[1]+.2), .005, (pos[3]-.5)]

    else:

        cabar_position = [(pos[2]+.005 +(0.055*(i%2))), (pos[1]-.22), .005, (pos[3]-.45)]

    cbar_ax = fig.add_axes(cabar_position)

    plt.colorbar(dense_map, cax=cbar_ax, orientation="vertical")

    # plt.colorbar(dense_map)

# fig.tight_layout()

y_hist.invert_xaxis()

x_hist.invert_yaxis()

img_name='multiple-dens_hist-map.png'

plt.savefig(img_name, dpi=300, bbox_inches='tight')

plt.show()
