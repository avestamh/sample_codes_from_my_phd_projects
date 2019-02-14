import os

def gen_sys_config(total_cores):
    for main_core in range(1, total_cores+1):
        for sys1_core in range(2, total_cores-main_core+1):
            for sys2_core in range(1, total_cores-main_core-sys1_core+1):
                if (main_core+sys1_core+sys2_core==total_cores) and (main_core<sys1_core) and (main_core<sys2_core) and (sys1_core>=sys2_core):
                    yield (main_core, sys1_core, sys2_core)


for gsc in  gen_sys_config(16):
    ntraj=3
    cmd_header = "sbatch --export "
    sys_config = "main_core=%d,sys1_core=%d,sys2_core=%d,"%(gsc[0], gsc[1], gsc[2])
    other_confg = "traj=%d,force=%d,machin=%s"%(5, 300, "clpy")
    cmd_str = cmd_header + sys_config + other_confg + " job.slurm"
    print cmd_str
    # os.system(cmd_str)
