#!/bin/bash

###############################
#### 256x192
###############################

#sed -i -e "14s/!integer/ integer/g" digda_mpi.f90
#sed -i -e "16s/ integer/!integer/g" digda_mpi.f90

#sed -i -e "100,102s/!integer/ integer/g" digda_grid.f90
#sed -i -e "110,112s/ integer/!integer/g" digda_grid.f90

#sed -i -e "s/digda_21601080/digda_256192/g" Mkinclude

 ##############################
 ### 2160x1080
 ##############################

 sed -i -e "14s/ integer/!integer/g" digda_mpi.f90
 sed -i -e "16s/!integer/ integer/g" digda_mpi.f90

 sed -i -e "100,102s/ integer/!integer/g" digda_grid.f90
 sed -i -e "110,112s/!integer/ integer/g" digda_grid.f90

 sed -i -e "s/digda_256192/digda_21601080/g" Mkinclude

 ##############################
 module purge
 module load NECSDK/all
 module load zlib_ve/1.2.11
 module load szip_ve/2.1.1
 module load NetCDF4_ve/netcdf-c/4.7.2
 module load NetCDF4_ve/netcdf-f/4.5.2
 module load wgrib2_netcdf4/3.0.2

 # compile
 sed -i -e "s/ES4s_intel/ES4v/g" Mkinclude
 make clean
 make > make.log_nfort 2>&1 &
 MAKE_PID=$!
 tail -f make.log_nfort &
 TAIL_PID=$!
 wait $MAKE_PID
 kill $TAIL_PID
 ###sed -i -e "s/ES4v/ES4s_intel/g" Mkinclude

