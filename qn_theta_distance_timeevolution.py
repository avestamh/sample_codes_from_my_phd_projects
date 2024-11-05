
# Written by: Sadra Avestan
# Date: November 13, 2019
import matplotlib.pyplot as plt

import numpy as np

import math

from sadra import pubfig as pf

pf.setup(28, 18)

pf.setup(label_font=30, tick_font=28, axis_width=2, tick_major_width=2, tick_minor_width=1, tick_major_size=4, tick_minor_size=2, showminorticks=True)

force_dir = "address1"

angle_dir = "address2"

qn_dir = "address3"

smooth_scale = 10

smooth_window = np.ones(smooth_scale)/smooth_scale

param_name = ['force', 'QN', 'theta', 'PHI']

CrciticalForce = []

for traj in range(1, 105):

    qn_address =  qn_dir + "traj%s/qn-gfp-strands-k_5-rate_0.005-t300.dat"%traj

    force_file = force_dir + "/pulling-go-6atpase-s3-gfp-rate0.005-k5-300-traj%s/smd5.2204-2454_1.juj"%traj

    angle_name = angle_dir + "/1gfp-inertia-c_pull-rate0.005-k5-T300-l22-traj%s/1gfp-c_pull-cv-mom_inertia-.dat"%traj

    # get the appropriate length

    qn = np.genfromtxt(qn_address)[:,1]

    angle = np.genfromtxt(angle_name)

    theta = np.genfromtxt(angle_name)[:,1]

    force = np.genfromtxt(force_file)[:, 3]

    arr_len = min(len(force), len(theta), len(qn))

    # print arr_len

    ### QN

    qn = np.genfromtxt(qn_address)[:arr_len,1]

    QN = np.convolve(qn, smooth_window, mode='valid')

    # THETA

    theta = angle[:arr_len,1]

    theta_len = int(len(theta))

    THETA = np.convolve(theta, smooth_window, mode='valid')

    ## CALCULATE DISTANCE OF C-MASS FROM Z AXIS

    Xcor = angle[:,5]

    Ycor = angle[:,6]

    Zcor = angle[:,7]

    dofcofMass_zax = np.sqrt(Xcor**2 + Ycor**2)[:arr_len]

    DofCofMass_Zax = np.convolve(dofcofMass_zax, smooth_window, mode='valid')

    # PHI

    X = angle[:,2]  # principal moment of inertia

    Y = angle[:,3]

    # Z = angle[:,4]*0 if I want just to calclulate after projecting to the x-y surface

    Z = angle[:,4]

    # print Z

    I = np.array([X, Y, Z])

    # print I.shape

                ## *** to distinguish b/w two sp orientations whith 180 degrees difference in configuratrion but gives similar PHI ***

    Xnter = angle[:, 8]

    Ynter = angle[:, 9]

    Znter = angle[:, 10]

                # dfine a vevtor for comparison

    V = np.array([(Xcor-Xnter), (Ycor-Ynter), (Zcor-Znter)])

    # print V.shape

    V_dot_I = np.multiply(V, I).sum(axis=0)

    # print type(V_dot_I)

    # print V_dot_I.shape

    coef = (-1)**(V_dot_I < 0) ## to return 1 and -1 

    # for v, i in zip(V_dot_I, coef):

    #     print v, i

    X = np.multiply(X, coef)  ## element-wise multiplication

    Y = np.multiply(Y, coef)

                #****************************************

    PHI = (np.arccos(X/np.sqrt(X**2+Y**2))*180/np.pi + (Y<0)*180)[:arr_len]

    PHI = np.convolve(PHI, smooth_window, mode='valid')

    # FORCE

    force = np.genfromtxt(force_file)[:arr_len,3]

    FORCE = np.convolve(force, smooth_window, mode='valid') 

    ### MAX FORCE

    CrciticalForce.append(force.max())

    # TIME

    time_fact = 0.1

    TIME = [(i+1)*time_fact * 5 for i in range(len(PHI))] 

    # Plot

    fig, axs = plt.subplots(5, 1)

    axs[0].plot(TIME, FORCE, color='tab:red')

    axs[0].set_ylim(-30, 160)

    axs[0].set_ylabel('F (pN)')

    axs[1].plot(TIME, QN, color='tab:blue')

    axs[1].set_ylim(-0.2, 1.2)

    axs[1].set_ylabel('$Q_N$')

    axs[2].plot(TIME, THETA, color='tab:green')

    axs[2].set_ylim(2, 92)

    axs[2].set_ylabel(r'$\theta$ ($^\circ$)')

    axs[3].plot(TIME, PHI, color='tab:olive')

    axs[3].set_ylim(-20, 360)

    axs[3].set_xlabel(r't ($\tau$)')

    axs[3].set_ylabel(r'$\phi$ ($^\circ$)')

    axs[4].plot(TIME, DofCofMass_Zax, color='tab:pink')

    axs[4].set_ylim(0, 35)

    axs[4].set_xlabel(r't ($\tau$)')

    axs[4].set_ylabel(r'd ($\AA$)')

    # Adjust

    plt.subplots_adjust(hspace=0)

    for i in range(5):

        # if i < 3:

        #     pass

            # axs[i].set_xticks([])

        axs[i].set_xlim(min(TIME), max(TIME))

        # axs[i].grid(True, which='both', linestyle='--')

        axs[i].grid(True, which='major', c='dimgray', linestyle='-')

        axs[i].grid(True, which='minor', linestyle='--')

    for ax in axs:

        ax.yaxis.set_label_coords(-.11, 0.5)  ## align the y axes

    fig.set_size_inches(10,10)

    plt.savefig("img-time-evo-traj%d.png"%traj, bbox_inches='tight', dpi=300)

    plt.close()

### SAVE THE CRITICAL FORCE

# with open('CrciticalForce.dat', 'w+') as file:

#     for listitems in CrciticalForce:

#         file.write("%s\n"%listitems)

# print len(CrciticalForce)
