#!/bin/bash
#SBATCH --ntasks=150
#SBATCH --nodes=5
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3GB
#SBATCH --time=1:00:00
#SBATCH --account=jpwalter_533
#SBATCH --partition=epyc-64

module load gcc/8.3.0
module load openmpi/4.0.2
module load pmix

mpicc -o Lab5 Lab5.c
export OMP_NUM_THREADS=20
mpirun -np 150 ./Lab5
