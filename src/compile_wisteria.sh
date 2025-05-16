#!/bin/bash

#resolution="128x64"
#resolution="256x192"
 resolution="2160x1080"

 ##############################
 ### 128x64
 ##############################

 if [ "$resolution" = "128x64" ]; then

 sed -i -e "10,23s/ integer/!integer/g" digda_mpi.f90
 sed -i -e "10s/!integer/ integer/g"    digda_mpi.f90 ; ncpu="8"

 sed -i -e "164,166s/!integer/ integer/g" digda_grid.f90
 sed -i -e "169,171s/ integer/!integer/g" digda_grid.f90
 sed -i -e "179,181s/ integer/!integer/g" digda_grid.f90

 sed -i -e "s/digda_256192/digda_128064/g"   Mkinclude
 sed -i -e "s/digda_21601080/digda_128064/g" Mkinclude

 fi

###############################
#### 256x192
###############################

 if [ "$resolution" = "256x192" ]; then

 sed -i -e "10,23s/ integer/!integer/g" digda_mpi.f90
 sed -i -e "16s/!integer/ integer/g"    digda_mpi.f90 ; ncpu="48"

 sed -i -e "164,166s/ integer/!integer/g" digda_grid.f90
 sed -i -e "169,171s/!integer/ integer/g" digda_grid.f90
 sed -i -e "179,181s/ integer/!integer/g" digda_grid.f90

 sed -i -e "s/digda_128064/digda_256192/g"   Mkinclude
 sed -i -e "s/digda_21601080/digda_256192/g" Mkinclude

 fi

###############################
#### 2160x1080
###############################

 if [ "$resolution" = "2160x1080" ]; then

 sed -i -e "10,23s/ integer/!integer/g" digda_mpi.f90
 sed -i -e "18s/!integer/ integer/g"    digda_mpi.f90 ; ncpu="90"

 sed -i -e "164,166s/ integer/!integer/g" digda_grid.f90
 sed -i -e "169,171s/ integer/!integer/g" digda_grid.f90
 sed -i -e "179,181s/!integer/ integer/g" digda_grid.f90

 sed -i -e "s/digda_128064/digda_21601080/g" Mkinclude
 sed -i -e "s/digda_256192/digda_21601080/g" Mkinclude

 fi
 ##############################

 sed -i -e "9s/digda/digda_${resolution}.${ncpu}/g" Mkinclude

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

 ##############################

 sed -i -e "9s/digda_${resolution}.${ncpu}/digda/g" Mkinclude

 ##############################
