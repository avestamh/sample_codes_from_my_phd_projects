
#-------------------------------------------------------------------------------
# Written by: Sadra Avestan
# Date: April 22, 2020
# Description: calculating of the average hydrogen bonds in all subunits of ClpP
#---------------------------------------------------------------------------------
import numpy as np

import matplotlib.pyplot as plt

from sadra import pubfig as pf

pf.setup()

ntraj = 4

avrg = np.zeros((500,14))

for traj in range(1,ntraj+1):

    all_data = np.genfromtxt('clpp_close-hbond-traj%d/hbond.dat'%traj, comments='#')[:,1:15]

    avrg = avrg+all_data

avrg=avrg/4

# np.savetxt('hbond-avrg-allloop.txt', avrg, fmt='%s')

fig, axs = plt.subplots(2,1)

cis_data = avrg[:,0:7].T

trans_data = avrg[:,7:14].T

im = axs[0].imshow(cis_data, interpolation='nearest', cmap='jet', vmin=0, vmax=3)

im = axs[1].imshow(trans_data, interpolation='nearest', cmap='jet', vmin=0, vmax=3)

plt.subplots_adjust(hspace=0.3)

for i,ax in enumerate(axs):

    ax.set_xlim(-.5, 49.5)

    xticks = [(i+1) for i in range(50)]

    # xlabs = [(i+1) for i in range(0,50,100)]

    ax.set_yticks([0, 1, 2, 3, 4, 5, 6])

    if i>0:

        ax.set_xlabel('Time (ns)')

        ax.xaxis.labelpad = -3

    if i>0:

        ax.set_yticklabels(['H', 'I', 'J', 'K', 'L', 'M', 'N'], ha='center')

    else:

        ax.set_yticklabels(['A', 'B', 'C', 'D', 'E', 'F', 'G'], ha='center')    

    for j in range(8):

        ax.axhline(y=j-0.5, color='black')

    xmin, xmax = ax.get_xlim()   

    ax.set_aspect(xmax*5e-2)

    ax.tick_params(axis='y', which='major', pad=10)

    ax.tick_params(axis='y', which='minor', left=False)

    ax.tick_params(axis='y', which='minor', right=False)

cbar = fig.colorbar(im, ax=axs.ravel().tolist(), shrink=1.0, aspect=30)

cbar.set_ticks([1, 2, 3])

import string

subfig_labels = string.ascii_uppercase

for i, ax in enumerate(axs.flatten()):

     ax.text(0.00, 1.05, "(%s)"%subfig_labels[i], transform=ax.transAxes, size=16)

pf.save3('habond-avrg-alltraj')
