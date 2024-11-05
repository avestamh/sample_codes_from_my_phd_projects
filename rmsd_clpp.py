#---------------------------------------------------------------
# Written by: Sadra Avestan
# Date: March 7 2021
# Description: calculation of rmsd for ClpP
#-----------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

from sadra import pubfig as pf
pf.setup(label_font=20, tick_font=20, axis_width=2, tick_major_width=3,
         tick_minor_width=2.5, tick_major_size=4, tick_minor_size=3, showminorticks=True)

plt.rcParams['figure.figsize'] = (5, 10)
file_address = 'clpp_close-rmsd'
mutations = ['wild-type', 'e14a-r15a', 'e8k', 'k25e',
        'e14a-r15a-e8a', 'e14a-r15a-k25a', 'i7p']


ntraj = 4
time_fact = 0.1
time = [(i+1)*time_fact for i in range(500)]

fig, axs = plt.subplots(7,1)
fig.subplots_adjust(hspace=.3)


for i,mute in enumerate(mutations):
    rmsd_average = np.zeros(500,)
    for traj in range(1,ntraj+1):
        rmsd_file = file_address +'/clpp_close-rmsd-%s-traj%d/rmsd.dat'%(mute, traj)
        rmsd_data = np.genfromtxt(rmsd_file, comments="#")[:,1]


        axs[i].plot(time, rmsd_data, c='gray', alpha=.3) ### plot for the single trajs

        rmsd_average = rmsd_average + rmsd_data
        print('for traj: {} of mute :{}'.format(traj, mute))
        print(rmsd_data[100:103])
    rmsd_average = rmsd_average/ntraj

    print('the average for mute: {}'.format(mute))

    print(rmsd_average[100:103])
    print('*'*30)
    axs[i].plot(time, rmsd_average, c='r', lw=2)
    axs[i].set_title(mutations[i])
    axs[i].set_ylim(1,3)
    axs[i].set_xlim(0,50)
    axs[3].set_ylabel(r'RMSD ($\AA$)')
    if i<6:
        axs[i].set_xticklabels([])
    else:
        axs[i].set_xlabel('Time (ns)')

        # axs[i].ticK_params(axis='x', which='both', top=False)

plt.savefig('rmsd.png', bbox_inches='tight', dpi=300)
# plt.show()
