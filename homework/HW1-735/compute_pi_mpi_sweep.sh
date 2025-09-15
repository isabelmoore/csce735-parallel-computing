#!/bin/bash
#SBATCH -J pi_mpi_param
#SBATCH -o out.%j
#SBATCH -e out.%j
#SBATCH -p short
#SBATCH -t 00:15:00
#SBATCH --mem=8G

set -euo pipefail
module load intel

# compile (safe to compile inside job; or comment out if already built)
mpiicx -O3 -xHost -o compute_pi_mpi.exe compute_pi_mpi.c

N=10000000000     # n = 1e10
P=64              # total MPI ranks

echo "tpn=${SLURM_NTASKS_PER_NODE}  nodes=${SLURM_JOB_NUM_NODES}  ntasks=$P"
# Run and capture exactly the time number your program prints.
# Your program prints: ... time (sec) =   2.2093
LINE=$(srun -n "$P" ./compute_pi_mpi.exe "$N" | awk -F'=' '/time \(sec\)/{gsub(/ /,"",$2); print $2}')

echo "RESULT tpn=${SLURM_NTASKS_PER_NODE} time_sec=${LINE}"

