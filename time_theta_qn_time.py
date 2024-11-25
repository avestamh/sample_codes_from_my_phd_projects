#!/usr/bin/env python
import bio.physics
# import matplotlib as mpl
# print mpl.__path__
import matplotlib.pyplot as plt
import numpy as np
bio.physics.sucks(14)

N = 100
y_axis = ["qn", "theta"]
for each_y in y_axis:
# fig, ax = plt.subplots(2, 5)  
    for traj in range(6, 9):
        # if traj==5:
        #     continue
        qn_name = "qn-traj%d/qn.dat.merged"%traj
        theta_name = "theta-traj%d/theta.dat.merged"%traj
        
        qn_file = np.genfromtxt(qn_name)
        theta_file = np.genfromtxt(theta_name)
        qn_data = qn_file[:,3]
        theta_data = theta_file[:,3]
        # time = len(qn_data)
        # timefact = time*0.01 # 
        if each_y == "theta":
            y = theta_data
            plt.xlim(0, 90)
            plt.ylabel("$\\theta$ [$^\circ$]", fontsize=20)
        else:
            y = qn_data
            plt.ylim(0, 1)
            plt.ylabel("$Q_N$", fontsize=16)

        y = np.convolve(y, np.ones((N,))/N, mode="valid")
        time = np.cumsum(np.ones(y.shape)*0.01)
        plt.plot(time, y)
        plt.xlabel("time (ns)", fontsize=16)  
        img_name = "scatter-plot-%s-time-traj%d.png"%(each_y, traj)
        plt.savefig(img_name, format='png', dpi=300, bbox_inches='tight')
        plt.close()
