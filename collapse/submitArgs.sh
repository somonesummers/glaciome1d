#!/bin/bash 
#SBATCH -J L6000_slow # job name 
#SBATCH -o output_%j.txt # output and error file name (%j expands to jobID)
#SBATCH --account=gts-arobel3-atlas    #charge account
#SBATCH -N1 --ntasks-per-node=1   #total number of nodes,CPUs requested
#SBATCH --mem-per-cpu=1G
#SBATCH -qinferno
#SBATCH -t 60:10:49 # run time (hh:mm:ss)
#SBATCH --mail-user=psummers8@gatech.edu
#SBATCH --mail-type=end,fail  # email me when the job finishes/fails

cd $SLURM_SUBMIT_DIR    # Change to working directory
set -e
#module load python/3.10.10
#module spider anaconda3/2023.03
# bash ../makeRunMpi.sh
~/.conda/envs/MITgcm/bin/python -u ../runModelSweepEqualib.py -s 1 -u 6000 -m .35 -p .7
~/.conda/envs/MITgcm/bin/python -u ../runModelUnstable.py -s 1
