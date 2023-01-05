#!/bin/bash
#SBATCH --constraint=intel
#SBATCH   --partition=contrib            # submit   to the normal(default) partition
#SBATCH   --job-name=TaC-test             # name the job
#SBATCH   --output=TaC.out        # write stdout/stderr   to named file
#SBATCH   --error=TaC-%j.err      
##SBATCH   --time=0-02:00:00             # Run for max of 02 hrs, 00 mins, 00 secs
#SBATCH   --nodes=2
#SBATCH --ntasks-per-node=48                    # Request N nodes
##SBATCH   --cpus-per-task=16            # Request n   cores per node
##SBATCH   --mem-per-cpu=2GB             # Request nGB RAM per core
#SBATCH --mem=10GB
#SBATCH --export=ALL

#load modules with  
module load python
module load gnu10
module load OneAPI/2022.1.2  
mpirun -np 96 /home/ikhatri/local/bin/vasp.std.hopper >stdout
