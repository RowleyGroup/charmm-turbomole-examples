#!/usr/bin/env bash
#PBS -q qtest
#PBS -l walltime=0:10:00
#PBS -l nodes=1:ppn=1
#PBS -N indole_ex
#PBS -o std.out
#PBS -j oe

#load the compilers for CHARMM
module load gcc/4.8.2
cd $PBS_O_WORKDIR

/home/bin/c38b1/exec/gnu/charmm < charmm_input/indole_ex.inp > charmm_output/indole_ex.out
