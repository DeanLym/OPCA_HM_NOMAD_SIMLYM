#!/bin/bash
#PBS -N OPCA_HM_NOMAD_SIMLYM
#PBS -l nodes=9:ppn=16
#PBS -q default
#PBS -V
#PBS -m e
cd $PBS_O_WORKDIR
#
mpirun -np 144 OPCA_HM_NOMAD_SIMLYM  >> out.txt

# end script
