import numpy as np 

import matplotlib.pyplot as plt 

import scipy

import math

from mpl_toolkits.mplot3d import Axes3D

from scipy.stats import gaussian_kde

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from operator import itemgetter as ig

from pylab import rcParams

from sadra import pubfig as pf

# pf.setup(18,16)

pf.setup(label_font=28, tick_font=20, axis_width=2, tick_major_width=2, tick_minor_width=2, tick_major_size=4, tick_minor_size=3, showminorticks=True)

angle_path = '/media/data0/sadra/26s-S3/atpaseparts/gomodel/2gfp/analysis/editedparam/2gfp-rep-force/inertia-azimuth'

qn_path    = '/media/data0/sadra/26s-S3/atpaseparts/gomodel/2gfp/analysis/editedparam/2gfp-rep-force/repetitve_random-qn-rmsd-rg-2gfp-26s-l51-T300-f100/repetitive-rand_pull-6atpase-nocyc400'

patha = [1, 2, 4, 5, 6, 8, 10, 12, 14, 15, 16, 18, 20, 23, 26, 28, 29, 31, 33, 34, 36, 38, 39, 40, 43, 44, 45, 46, 48, 49, 51, 55, 56, 57, 62, 63, 64, 66]

pathb = [3, 7, 9, 17, 19, 21, 22, 25, 27, 30, 32, 37, 41, 42, 47, 50, 52, 53, 54, 58, 59, 61, 65]

path_name = ['patha', 'pathb']

prot_name = ['ngfp', 'gfpc']

prot_column = {'ngfp':3, 'gfpc':4}

def get_theta_phi(path, trajs):

        ## phi

    P_ngfp = []

    P_gfpc = []

    ## theta

    T_ngfp = []

    T_gfpc = []

    for traj in trajs: 

        for prot in prot_name:

            theta_file = angle_path+'/%s-rept-inertia-T300-l51-f100/%s-repet-first_inertia-f100-traj%d/%s-rep-f100-mom_inertia.dat.merged'%(prot,prot, traj, prot)

            qn_file = qn_path + "/qn-rmsd-2gfp-26s-repetitive-l51-f100-traj%d/qn-rept-2gfp.dat.merged"%traj

            qn_data = np.genfromtxt(qn_file)[:,prot_column[prot]]

            qn_8 = np.where(qn_data<0.8)[0][0]

            # print(qn_8)

            inp_theta = np.genfromtxt(theta_file)

            theta_data = inp_theta[0:qn_8,3]

                        # principal moment of inertia

            X = inp_theta[0:qn_8,4]

            Y = inp_theta[0:qn_8,5]

            Z = inp_theta[0:qn_8,6]

            ## coors of COM

            Xcor = inp_theta[0:qn_8,7]

            Ycor = inp_theta[0:qn_8,8]

            Zcor = inp_theta[0:qn_8,9]

                          ## *** to distinguish b/w two sp orientations whith 180 degrees difference in configuratrion but gives similar PHI ***

            Xnter = inp_theta[0:qn_8, 10]

            Ynter = inp_theta[0:qn_8, 11]

            Znter = inp_theta[0:qn_8, 12]

            I = np.array([X, Y, Z]) # vector of moment of inertia

                            # dfine a vevtor for comparison

            V = np.array([(Xcor-Xnter), (Ycor-Ynter), (Zcor-Znter)])

            V_dot_I = np.multiply(V, I).sum(axis=0)

            coef = (-1)**(V_dot_I<0)

            #     for v, i in zip(V_dot_I, coef):

            #             print v, i

            X = np.multiply(X, coef)

            Y = np.multiply(Y, coef)

            PHI = np.arccos(X/np.sqrt(X**2+Y**2))*180/np.pi + (Y<0)*180 

            T_ngfp.extend(theta_data) if prot == 'ngfp' else T_gfpc.extend(theta_data)

            P_ngfp.extend(PHI) if prot == 'ngfp' else P_gfpc.extend(PHI)

        ##__________________New X-axis w/r to the coiled coil region of SEGM and SEGDL of proteasome___________________________

            ### extracting new angle

    NewXvec = np.array([35.4941115, 32.317318 , 0])

    mainX = np.array([1, 0, 0])

    Vec = np.dot(NewXvec, mainX)

    NewXvecSize = np.sqrt((NewXvec[0])**2 + (NewXvec[1])**2) ## The size of the x vec is 1, not included

    CCPhi = np.arccos(Vec/NewXvecSize) * 180/np.pi

    Tn = np.array(T_ngfp)

    Tc = np.array(T_gfpc)

    Pn = np.array(P_ngfp)

    Pc = np.array(P_gfpc)

    Pn = Pn - CCPhi

    Pn = Pn + (Pn<0)*360

    Pc = Pc - CCPhi

    Pc = Pc + (Pc<0)*360

    return Tn, Pn, Tc, Pc

Tna, Pna, Tca, Pca = get_theta_phi('patha', patha)

Tnb, Pnb, Tcb, Pcb = get_theta_phi('pathb', pathb)

# print(len(Tna),Tna,len(Tnb), len(Pnb), Tnb)

import scipy.stats as st

### _____generate the meshgrid plots__________

# dummy functions to setup the limit

def gen_lowerboundary(num):

    return num*1.01 if num < 0 else num*0.99

def gen_upperboundary(num):

    return num*0.99 if num < 0 else num*1.01 

def Kdeplot(axscolln, idx, fig, x, y, nlevels=10, nbins=100):

    ## ----------Generate limits for plot and meshgrid---------

    xmin, xmax = gen_lowerboundary(x.min()), gen_upperboundary(x.max())

    ymin, ymax = gen_lowerboundary(y.min()), gen_upperboundary(y.max())

    ### ----run meshgrid and  run kernel density stimation-------Evaluate kde on a grid

    xMesh, yMesh = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]

    positions = np.vstack([xMesh.ravel(), yMesh.ravel()]) #                                  #positions = np.vstack([item.ravel() for item in [xMesh, yMesh]])

    values = np.vstack([x, y])

    kernel = st.gaussian_kde(values)

    zMesh = np.reshape(kernel(positions).T, xMesh.shape)

    levels = axscolln[idx].contour(xMesh, yMesh, zMesh, nlevels, colors='k', alpha=0).levels

    cfplot = axscolln[idx].contourf(xMesh, yMesh, zMesh, levels[1:], cmap='jet')

    axins = inset_axes(axscolln[idx], width="4%", height="30%", loc='upper left', borderpad=.9)

    cbar = fig.colorbar(cfplot, cax=axins, orientation="vertical",format="%1.1e", ticks=ig(1, -1)(levels))

    cbar.ax.tick_params(labelsize=12) 

    cbar.ax.minorticks_off()

    if idx == 0:

        axscolln[idx].set_ylabel(r'$\phi$ ($^\circ$)')

        axscolln[idx].set_xlabel(r'$\theta$ ($^\circ$)')

    else:

        axscolln[idx].set_xlabel(r'$\theta$ ($^\circ$)')

    axscolln[idx].set_xlim(0, 90)

    axscolln[idx].set_ylim(0, 360)

    name = "Major-NGFP" if idx==0 else "Major-GFPC" if idx==1 else "Minor-NGFP" if idx==2 else "Minor-GFPC"

    axscolln[idx].set_title(name)

    axscolln[idx].set_xticks([i*25 for i in [0, 1, 2, 3, 3.6]])

    axscolln[idx].set_xticklabels([0, 25, 50, 75])

    if idx==0:

        axscolln[idx].tick_params(axis='y', which='both', right=False)

    if idx==1 or idx==2:

        axscolln[idx].set_yticks([])

    if idx==3:

        axscolln[idx].set_yticklabels([])

        axscolln[idx].tick_params(axis='y', which='both', left=False)    

rcParams['figure.figsize']=8,5

fig, axs = plt.subplots(1,4)

plt.subplots_adjust(wspace=0)

Kdeplot(axs, 0, fig, Tna, Pna)

Kdeplot(axs,1,fig, Tca, Pca)

Kdeplot(axs,2,fig, Tnb, Pnb)

Kdeplot(axs,3,fig, Tcb, Pcb)

# plt.show()

img_name = 'rrf-theta-phi.png'

plt.savefig(img_name, dpi=300, bbox_inches='tight')
