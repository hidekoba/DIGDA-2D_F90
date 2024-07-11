#!/bin/bash

#module avail
 module unload intel
 module unload impi
 module load fj
 module load fjmpi
 module load hdf5
 module load phdf5
 module load netcdf
 module load netcdf-fortran

 make clean
 make

 module unload fj
 module unload fjmpi
 module load intel
 module load impi
