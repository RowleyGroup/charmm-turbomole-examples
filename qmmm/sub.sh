#!/usr/bin/env bash
#PBS -q qtest
#PBS -l walltime=0:10:00
#PBS -l nodes=1:ppn=1
#PBS -N mg(h2o)6
#PBS -o std.out
#PBS -j oe

#load the compilers for CHARMM
module load gcc/4.8.2
cd $PBS_O_WORKDIR

/home/bin/c38b1/exec/gnu/charmm < charmm_input/mg_6h2o_qmmm.inp > charmm_output/mg_6h2o_qmmm.out
