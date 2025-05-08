#!/bin/bash

###############################
#### 256x192
###############################

#sed -i -e "14s/!integer/ integer/g" digda_mpi.f90
#sed -i -e "16s/ integer/!integer/g" digda_mpi.f90

#sed -i -e "133,135s/!integer/ integer/g" digda_grid.f90
#sed -i -e "143,145s/ integer/!integer/g" digda_grid.f90

#sed -i -e "s/digda_21601080/digda_256192/g" Mkinclude

 ##############################
 ### 2160x1080
 ##############################

 sed -i -e "14s/ integer/!integer/g" digda_mpi.f90
 sed -i -e "16s/!integer/ integer/g" digda_mpi.f90

 sed -i -e "133,135s/ integer/!integer/g" digda_grid.f90
 sed -i -e "143,145s/!integer/ integer/g" digda_grid.f90

 sed -i -e "s/digda_256192/digda_21601080/g" Mkinclude

 ##############################

#module avail

 module unload netcdf
 module unload netcdf-fortran
 module unload intel
 module unload impi

 module load fj
 module load fjmpi
 module load hdf5
 module load phdf5
 module load netcdf
 module load netcdf-fortran

 module list

 make clean
 make > make.log 2>&1 &
 MAKE_PID=$!
 tail -f make.log &
 TAIL_PID=$!
 wait $MAKE_PID
 kill $TAIL_PID

#module unload fj
#module unload fjmpi
#module load intel
#module load impi

#module list
