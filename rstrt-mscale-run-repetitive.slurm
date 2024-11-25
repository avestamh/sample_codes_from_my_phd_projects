#!/bin/bash
#SBATCH --job-name="charmm"
#SBATCH -A cin108
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 48:00:00
#SBATCH --export=ALL

#-------------------------------------------------------------------------
export LD_LIBRARY_PATH=/home/wangqi/shared/gcc-5.2.0/lib64:$LD_LIBRARY_PATH
export PATH=/home/wangqi/shared/gcc-5.2.0/bin:$PATH

export PATH=/home/sadra/shared/openmpi/bin:$PATH
export MPI_LIB=/home/sadra/shared/openmpi/lib
export MPI_INCLUDE=/home/sadra/shared/openmpi/include
export LD_LIBRARY_PATH=/home/sadra/shared/openmpi/lib:$LD_LIBRARY_PATH
#-------------------------------------------------------------------------

#-----------------------variables------------
ntraj=$ntraj
T=300
fmean=$fmean
sigma=$sigma
domain=$domain
machin=$machin
maxcycle=500
level=mscale #go
fch=0.01 ## consharm force
maincore=$maincore
sys1core=$sys1core
sys2core=$sys2core
sprotate=$((RANDOM%360+1))
iseed=$RANDOM
jobname=$level-$machin-$domain
#-----------------------------------------------

#-----------------directions--------------------
tmpdir=/scratch/$USER/$SLURM_JOB_ID
jobdir=$jobname-f$fmean-nbfx10-2lp-nogeo-traj$ntraj
#jobdir=$jobname-f$fmean-2lp-geo-traj$ntraj
#jobdir=bmrk-$jobname-f$fmean-2lp-nogeo-m$maincore-so$sys1core-st$sys2core-traj$ntraj


inpdir=${SLURM_SUBMIT_DIR}
psfdir=/home/sadra/i27-clpy/psfgen/toppar
outdir=/oasis/projects/nsf/cin108/sadra/mscal/i27-clpy/trajs
CHMEXE=/home/sadra/app/c40b1/exec/gnu_M/charmm.mpi.mscale.tmd.cvel.megatypes
#-----------------------------------------------
cd      ${tmpdir}
mkdir   ${jobdir}
cd      ${jobdir}


cp $CHMEXE 			charmm
cp $inpdir/*.inp		.
cp -r $psfdir			.

##_______________++RUN START AND THEN THE RESTAR FILE++_______________##

mpirun -np $maincore ./charmm icyc=$icyc maincore=$maincore sys1core=$sys1core sys2core=$sys2core  ntraj=$ntraj fmean=$fmean sigma=$sigma T=$T domain=$domain machin=$machin level=$level imode=$imode fch=$fch sprotate=$sprotate iseed=$iseed < $jobname-strt.inp > $jobname-strt.out
#wait


mode=(pull relax)

for icyc in `seq 1 $maxcycle`; do
	for imode in ${mode[@]}; do
		
		mv restart.res restart.rea
		
		mpirun -np $maincore ./charmm icyc=$icyc maincore=$maincore sys1core=$sys1core sys2core=$sys2core  ntraj=$ntraj fmean=$fmean sigma=$sigma T=$T domain=$domain machin=$machin imode=$imode level=$level fch=$fch < $jobname.inp > $jobname-$imode-$icyc.out
		wait
		gzip -f $jobname-$imode-$icyc.out
		if [ -e stop.flag ]; then break 2; fi
		#wait
	done
done
#----------------------------------------------------------------------

rm charmm # defined localy
rsync -auvz ${tmpdir}/${jobdir} ${outdir}/


###SBATCH --ntasks-per-node=24
