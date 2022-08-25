#!/bin/bash
#SBATCH --job-name=PVi90dgnP480t1e4yrs_P31yrs_fi0v995_fj0v995_dz_1v0kps
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --time=4320
#SBATCH --partition=long
#SBATCH --mail-user=mialy.rabenanahary@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=slurm-PVi90dgnP480t1e4yrs_P31yrs_fi0v995_fj0v995_dz_1v0kps-%A.out
#SBATCH --mem=64gb
#SBATCH --tmp=250gb

SCRATCH=/scratch/$USER/run.${SLURM_JOBID}.${SLURM_JOB_NAME}
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p $SCRATCH

#SCRATCH=/scratch/$USER/run

cd /data/$USER/Output
mkdir ${SLURM_JOB_NAME}

cd $SCRATCH
srun --ntasks=$SLURM_JOB_NUM_NODES cp /obs/mrabenanahary/Post-Traitement/geoCompiler_v2/test_input_i95.par .
srun --ntasks=$SLURM_JOB_NUM_NODES cp /obs/mrabenanahary/Post-Traitement/geoCompiler_v2/*.py .
srun --ntasks=$SLURM_JOB_NUM_NODES cp /obs/mrabenanahary/Post-Traitement/geoCompiler_v2/Jet_CC_Shang06_0296.vtu .
srun --ntasks=$SLURM_JOB_NUM_NODES cp /obs/mrabenanahary/Post-Traitement/geoCompiler_v2/Jet_CC_0100.vtu .
echo ${SLURM_JOBID}
echo $HOSTNAME
pwd

##srun --ntasks=$SLURM_JOB_NUM_NODES cp /data/$USER/MesDonnees .

#mpiexec ./MonProg > MonProg.out
#mpirun -n $SLURM_NTASKS amrvac -i zheng19.par
srun --ntasks=$SLURM_JOB_NUM_NODES more test_input_i90.par
#mpiexec amrvac -i ${SLURM_JOB_NAME}.par
srun --ntasks=$SLURM_JOB_NUM_NODES python3 geoCompiler.py ./test_input_i90.par
srun --ntasks=$SLURM_JOB_NUM_NODES mv ./${SLURM_JOB_NAME} /data/$USER/Output/${SLURM_JOB_NAME}/.
srun --ntasks=$SLURM_JOB_NUM_NODES cp /obs/mrabenanahary/Post-Traitement/geoCompiler_v2/slurm-${SLURM_JOB_NAME}-${SLURM_JOBID}.out /data/$USER/Output/${SLURM_JOB_NAME}/.
cd ${SLURM_SUBMIT_DIR}
#mv ${SCRATCH}/MonProg.out .
#srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf ${SCRATCH}

exit 0
