#!/bin/bash
#PBS -l ncpus=12
#PBS -l walltime=00:20:00
#PBS -l mem=3gb
#PBS -P ud25
#PBS -q normal

module purge
module load gcc openmpi
module load fftw3
cd $PBS_O_WORKDIR

##mpirun -np 2 /home/YOUR UNIQUE PATH TO YOUR EXE..... > output.log

##my job looks like this
mpirun -np 12 ./main  > ncpu_12_log.txt
