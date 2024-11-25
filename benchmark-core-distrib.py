import os

def gen_sys_config(totalcores):
    for maincore in range(2, totalcores+1):
        for sys1core in range(2, totalcores-maincore+1):
            for sys2core in range(2, totalcores-maincore-sys1core+1):
                if (maincore+sys1core+sys2core==totalcores) and (maincore<=sys1core) and (maincore<sys2core) and (sys1core<=sys2core):
                    yield (maincore, sys1core, sys2core)

def submit_jobs(totalcores):
    for gsc in  gen_sys_config(totalcores):
        cmd_header = "sbatch --ntasks-per-node=%d --export "%totalcores
        sys_config = "maincore=%d,sys1core=%d,sys2core=%d,"%(gsc[0], gsc[1], gsc[2])
        other_confg = "ntraj=%d,fmean=%d,sigma=%d,domain=%s,machin=%s"%(1, 200, 20, "i27", "1do2")
        cmd_str = cmd_header + sys_config + other_confg + " mscale-run-repetitive.slurm"
        print cmd_str
        os.system(cmd_str)
        # break

submit_jobs(8)
