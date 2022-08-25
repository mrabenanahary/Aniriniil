#!/bin/bash
#SBATCH --job-name=Lee2001_us_300yrs_0610_i90dgnP480
#SBATCH --nodes=1 --ntasks-per-node=48
#SBATCH --time=4320
#SBATCH --partition=long
#SBATCH --mail-user=mialy.rabenanahary@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=slurm-Lee2001_us_300yrs_0610_i90dgnP480-%A.out
#SBATCH --mem=64gb
#SBATCH --tmp=250gb

SCRATCH=/scratch/$USER/run.${SLURM_JOBID}.${SLURM_JOB_NAME}
srun --ntasks=$SLURM_JOB_NUM_NODES mkdir -p $SCRATCH

#SCRATCH=/scratch/$USER/run

cd /data/$USER/Output
mkdir ${SLURM_JOB_NAME}

cd $SCRATCH
srun --ntasks=$SLURM_JOB_NUM_NODES cp /obs/mrabenanahary/Post-Traitement/geoCompiler_v2/Lee2001_us_300yrs_0610yrs_i90_MV.par .
srun --ntasks=$SLURM_JOB_NUM_NODES cp /obs/mrabenanahary/Post-Traitement/geoCompiler_v2/*.py .
srun --ntasks=$SLURM_JOB_NUM_NODES cp /obs/mrabenanahary/Post-Traitement/geoCompiler_v2/Jet_CC_Shang06_0296.vtu .
srun --ntasks=$SLURM_JOB_NUM_NODES cp /obs/mrabenanahary/Post-Traitement/geoCompiler_v2/Lee2001_us_300yrs_0610.vtu .
echo ${SLURM_JOBID}
echo $HOSTNAME
pwd

##srun --ntasks=$SLURM_JOB_NUM_NODES cp /data/$USER/MesDonnees .

#mpiexec ./MonProg > MonProg.out
#mpirun -n $SLURM_NTASKS amrvac -i zheng19.par
srun --ntasks=$SLURM_JOB_NUM_NODES more Lee2001_us_300yrs_0610yrs_i90_MV.par
#mpiexec amrvac -i ${SLURM_JOB_NAME}.par
srun --ntasks=$SLURM_JOB_NUM_NODES python3 geoCompiler.py ./Lee2001_us_300yrs_0610yrs_i90_MV.par
srun --ntasks=$SLURM_JOB_NUM_NODES mv ./${SLURM_JOB_NAME} /data/$USER/Output/${SLURM_JOB_NAME}/.
srun --ntasks=$SLURM_JOB_NUM_NODES cp /obs/mrabenanahary/Post-Traitement/geoCompiler_v2/slurm-${SLURM_JOB_NAME}-${SLURM_JOBID}.out /data/$USER/Output/${SLURM_JOB_NAME}/.
cd ${SLURM_SUBMIT_DIR}
#mv ${SCRATCH}/MonProg.out .
#srun --ntasks=$SLURM_JOB_NUM_NODES rm -rf ${SCRATCH}

exit 0
