#!/bin/bash
#PBS -l walltime=00:10:00
#PBS -l select=2:ncpus=12:mpiprocs=12:mem=10000m
#PBS -m n
#PBS -N task_multMatrices

cd $PBS_O_WORKDIR
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')
echo "Number of MPI process: $MPI_NP"

mpicxx parallel.cpp -o parallel -std=c++11

echo "2 x 12 :"
mpirun -hostfile $PBS_NODEFILE -np $MPI_NP ./parallel 6168 3000 5184 2 12

echo "3 x 8 :"
mpirun -hostfile $PBS_NODEFILE -np $MPI_NP ./parallel 6168 3000 5184 3 8

echo "4 x 6 :"
mpirun -hostfile $PBS_NODEFILE -np $MPI_NP ./parallel 6168 3000 5184 4 6

echo "6 x 4 :"
mpirun -hostfile $PBS_NODEFILE -np $MPI_NP ./parallel 6168 3000 5184 6 4

echo "8 x 3 :"
mpirun -hostfile $PBS_NODEFILE -np $MPI_NP ./parallel 6168 3000 5184 8 3

echo "12 x 2 :"
mpirun -hostfile $PBS_NODEFILE -np $MPI_NP ./parallel 6168 3000 5184 12 2

