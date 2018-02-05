#!/bin/bash
# source me
module purge
module load intel/14.0.2.144
module load julia
module load R
module load openmpi/3.0.0
export PATH=/usr/local/julia/0.5.0_stand/julia-3c9d75391c/bin:$PATH
