#!/bin/bash
#PBS -l walltime=00:05:00
#PBS -l select=1:ncpus=8:mpiprocs=8:mem=10000m,place=scatter:exclhost
#PBS -m n
#PBS -N trace_2x4__analyzeMatrix

cd $PBS_O_WORKDIR

MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')
echo "Number of MPI process: $MPI_NP"

mpirun -trace -machinefile $PBS_NODEFILE -np $MPI_NP -perhost 2 ./parallel 5328 1300 4968 2 4

