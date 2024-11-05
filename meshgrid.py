import matplotlib.pyplot as plt

import numpy as np

from scipy.stats import gaussian_kde

import bio.physics

bio.physics.sucks(14)

# from sadra import test

# test.bark()

path_names = ['a', 'b']

prot_names = ['ngfp', 'gfpc']

beta_names = {'a':'b11', 'b':'b1'}

qn_columns = {'ngfp':3, 'gfpc':4}

def gen_lowerboundary(num):

    # return num

    return num*1.2 if num < 0 else num*0.8  # to set the boundary that we wnat to not be so compact (due to xmin and ymin)

def gen_upperboundary(num):

    # return num

    return num*0.8 if num < 0 else num*1.2

def kdeplot(x, y, nlevels=10, nbins=100):

    import numpy as np

    import matplotlib.pyplot as plt

    import scipy.stats as st

 #----- Generate limits for plot and meshgrid -----

    xmin, xmax = gen_lowerboundary(x.min()), gen_upperboundary(x.max())

    ymin, ymax = gen_lowerboundary(y.min()), gen_upperboundary(y.max())

    #----- Generate meshgrid and run kernel density estimation -----

    xMesh, yMesh = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]

    positions = np.vstack([xMesh.ravel(), yMesh.ravel()])

    values = np.vstack([x, y])

    kernel = st.gaussian_kde(values)

    zMesh  = np.reshape(kernel(positions).T, xMesh.shape)

    #----- Plot the results -----

    ax = plt.figure(figsize=(8,6)).gca()

    levels = ax.contour(xMesh, yMesh, zMesh, nlevels, colors='k', alpha=0).levels

    # levels = ax.contour(xMesh, yMesh, zMesh, nlevels, colors='k').levels

    levels = np.concatenate([[0.0],levels,[levels[-1]+levels[0]]])

    cfplot = ax.contourf(xMesh, yMesh, zMesh, levels[1:], cmap='jet') # to control the levels, that we want to keep or remove

    #----- Set up colorbar -----

    cbar = plt.colorbar(cfplot)

    cbar.ax.minorticks_off() # <----------------------------MOST IMPORTANT------------------------------>

    cbar.ax.xaxis.set_ticks_position('top')

    cbar.ax.xaxis.set_label_position('top')

    ticklabs = cbar.ax.get_yticklabels()

    cbar.ax.set_yticklabels(ticklabs,ha='right')

    cbar.ax.yaxis.set_tick_params(pad=40)

    cbar.set_label('Z', fontsize=16, labelpad=-45, y=1.05, rotation=0)

    # plt.title("2GFP CV %s"%prot)

    ax.set_xlabel("$\\theta$ ($^\circ$)", fontsize=20)

    ax.set_ylabel("$Q_N$", fontsize=20)

    ax.set_xlim(0, 90)

    ax.set_ylim(0, 1)

    # img_name = "path%s-prot-%s.png"%(path, prot)

    img_name = "new-mesggrid-plot-path%s-prot-%s.png"%(path, prot)

    plt.savefig(img_name, format='png', dpi=300, bbox_inches='tight', transparent=True)

    plt.close()

##############################

for path in path_names:

    for prot in prot_names:

        qn_name = "path%s-qn-%s-unified_with_theta_%s_0.2_rows.dat"%(path, prot, beta_names[path])

        theta_name = "path%s-%s-theta-0.2.dat"%(path, prot)

        print "Plotting for %s"%qn_name

        qn_file = np.genfromtxt(qn_name)

        theta_file = np.genfromtxt(theta_name)

        qn_data = qn_file[:,qn_columns[prot]]

        theta_data = theta_file[:,3]

        y  = qn_data

        x = theta_data

        kdeplot(x, y, nlevels=40, nbins=100)
