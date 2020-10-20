#!/bin/bash

#---------------------------------------------------------------------------------
# Account information

#SBATCH --account=staff              # basic (default), staff, phd, faculty

#---------------------------------------------------------------------------------
# Resources requested

#SBATCH --partition=gpu            # standard (default), long, gpu, mpi, highmem
#SBATCH --cpus-per-task=2          # number of CPUs requested (for parallel tasks)
#SBATCH --mem-per-cpu=4G           # requested memory
#SBATCH --time=0-8:00:00           # wall clock limit (d-hh:mm:ss)
#SBATCH --output=deepgmm.log       # join the output and error files

#---------------------------------------------------------------------------------
# Job specific name (helps organize and track progress of jobs)

#SBATCH --job-name=deepgmm         # user-defined job name

#---------------------------------------------------------------------------------
# Print some useful variables

echo "Job ID: $SLURM_JOB_ID"
echo "Job User: $SLURM_JOB_USER"
echo "Num Cores: $SLURM_JOB_CPUS_PER_NODE"

#---------------------------------------------------------------------------------
# Load necessary modules for the job

module load gcc/9.2.0
module load python/booth/3.6/3.6.3
module load mpi/openmpi-x86_64

#---------------------------------------------------------------------------------
# Commands to execute below...

ipython demand_experiment.py
