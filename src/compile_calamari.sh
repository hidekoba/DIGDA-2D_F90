#!/bin/bash

###############################
#### 256x192
###############################

#sed -i -e "14s/ integer/!integer/g" digda_mpi.f90
#sed -i -e "16s/ integer/1integer/g" digda_mpi.f90
#sed -i -e "10s/!integer/ integer/g" digda_mpi.f90

#sed -i -e "100,102s/!integer/ integer/g" digda_grid.f90
#sed -i -e "110,112s/ integer/!integer/g" digda_grid.f90

#sed -i -e "s/digda_21601080/digda_256192/g" Mkinclude

 ##############################
 ### 2160x1080
 ##############################

 sed -i -e "14s/ integer/!integer/g" digda_mpi.f90
 sed -i -e "16s/ integer/1integer/g" digda_mpi.f90
 sed -i -e "10s/!integer/ integer/g" digda_mpi.f90

 sed -i -e "100,102s/ integer/!integer/g" digda_grid.f90
 sed -i -e "110,112s/!integer/ integer/g" digda_grid.f90

 sed -i -e "s/digda_256192/digda_21601080/g" Mkinclude

 ##############################

 make clean
 make

