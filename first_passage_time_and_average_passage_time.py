## first passage time (FPT) and  average FPT

import numpy as np 
import matplotlib.pyplot as plt

patha = [1, 2, 4, 5, 6, 8, 10, 12, 14, 15, 16, 18, 20, 23, 26, 28, 29, 31, 33, 34, 36, 38, 39, 40, 43, 44, 45, 46, 48, 49, 51, 55, 56, 57, 62, 63, 64, 66]
pathb = [3, 7, 9, 17, 19, 21, 22, 25, 27, 30, 32, 37, 41, 42, 47, 50, 52, 53, 54, 58, 59, 61, 65]

path_name = ["a", "b"]
prot_name = ["ngfp", "gfpc"]
qn_columns = {"ngfp":3, "gfpc":4}

N = 20
overall_fpt = []

def write_array(inp_prot, inp_arr):
    avrg = np.mean(inp_arr)
    stddev = np.std(inp_arr)
    out_name = "fpt-stat-%s-avrg-%.2f-stddev-%.4f.dat"%(inp_prot, avrg, stddev)
    print out_name

    with open(out_name, 'w') as out_file:
        for ia in inp_arr:
            out_file.write("%.2f\n"%ia) 

    plt.hist(inp_arr, alpha=0.6, label=inp_prot, bins=10, color="#FF5722")
    plt.legend()
    plt.savefig('hist-%s.png'%inp_prot)
    plt.close()

def get_traj_fpt(trajs, path):
    fpts_ngfp = [] # first we should open an array to save everthing here
    fpts_gfpc = [] # if it were outside the def, it would accumulate for all the other and won't clear it each time the loop runs
    for traj in trajs:
       for prot in prot_name:
            inp_name = "repetitive-rand_pull-6atpase-nocyc400/qn-rmsd-2gfp-26s-repetitive-l51-f100-traj%s/qn-rept-2gfp.dat.merged"%traj
            inp_file = np.genfromtxt(inp_name)
            inp_data =  inp_file[:,qn_columns[prot]]
            time_fact = 0.01
            inp_data_smooth = np.convolve(inp_data, np.ones((N,))/N, mode='valid')
            # time = [(i+1)*time_fact for i in range(len(inp_data))]
            time = [(i+1)*time_fact for i in range(len(inp_data_smooth))]
            threshold=0.8
            fpt = 0.0 ## to make it global
            # for t, q in zip(time, inp_data):
            for t, q in zip(time, inp_data_smooth):
                if q < threshold:
                    fpt = t
                    break   ## so that it just print the first point not going through the loop

            if prot == 'ngfp':
                fpts_ngfp.append(fpt)
            else:
                 fpts_gfpc.append(fpt)

    write_array("ngfp-%s"%path, fpts_ngfp)
    write_array("gfpc-%s"%path, fpts_gfpc)

    return fpts_ngfp, fpts_gfpc

fn, fc = get_traj_fpt(patha, "patha")
overall_fpt.append(fn)
overall_fpt.append(fc)
fn, fc = get_traj_fpt(pathb, "pathb")
overall_fpt.append(fn)
overall_fpt.append(fc)
xlabel = "FPT"
xticks = [i*10 for i in range(8)]
labels = ['ngfp patha', 'gfpc patha', 'ngfp pathb', 'gfpc pathb']
colors = ['b', 'r', 'y', 'g']
fig, axs = plt.subplots(2,2)

for i, each_fpt in enumerate(overall_fpt):
    x = i/2
    y = i%2
    axs[x, y].hist(each_fpt, color=colors[i], label=labels[i])
    axs[x, y].legend()
    axs[x, y].set_xticks(xticks)
    axs[x, y].set_ylim(0, 10)
    axs[x, y].set_xlabel(xlabel)
    # plt.xlim(min(), max())

# plt.xticks(fontsize = 14)
# plt.yticks(fontsize = 14)
plt.xticks(np.arange(0, max(each_fpt), 10))

img_name = 'unfolding_time_66trajs.png'

plt.savefig(img_name, format='png', dpi=300, bbox_inches='tight')
