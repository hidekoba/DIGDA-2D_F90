#!/bin/bash
 set -euo pipefail

  resolution="128x64"
 #resolution="256x192"
 #resolution="2160x1080"

 #ciso=""
  ciso="_ciso"

 # ---------------------------------------------------------
 # Utility functions
 # ---------------------------------------------------------

 # Return line number of the first match of a regex in a file
 anchor_line () {
   local file="$1" regex="$2"
   grep -n -m1 -E "$regex" "$file" | cut -d: -f1
 }

 # Apply sed to a range relative to an anchor
 # sed_range <file> <regex> <offset_start> <offset_end> <sed-subst>
 sed_range () {
   local file="$1" regex="$2" off_s="$3" off_e="$4" subst="$5"
   local a s e
   a=$(anchor_line "$file" "$regex")
   s=$((a + off_s))
   e=$((a + off_e))
   echo "sed -i -e \"${s},${e}${subst}\" \"$file\""
   sed -i -e "${s},${e}${subst}" "$file"
 }

 # Apply sed to a single line relative to an anchor
 # sed_line <file> <regex> <offset> <sed-subst>
 sed_line () {
   local file="$1" regex="$2" off="$3" subst="$4"
   local a l
   a=$(anchor_line "$file" "$regex")
   l=$((a + off))
   echo "sed -i -e \"${l}${subst}\" \"$file\""
   sed -i -e "${l}${subst}" "$file"
 }

 # ---------------------------------------------------------
 # Anchors (patterns to detect base lines)
 # ---------------------------------------------------------
 ANCH_MPI='integer, parameter :: ncpu'
 ANCH_GRID='integer, parameter :: nxdim'
 ANCH_MK='bin/digda'

 # ---------------------------------------------------------
 # Replacement logic for each resolution
 # ---------------------------------------------------------
 case "$resolution" in
   "128x64")
     # digda_mpi.f90
     sed_range digda_mpi.f90 "$ANCH_MPI" 0 13 's/ integer/!integer/g'
     sed_line  digda_mpi.f90 "$ANCH_MPI" 0    's/!integer/ integer/g'
     ncpu="8"

     # digda_grid.f90
     sed_range digda_grid.f90 "$ANCH_GRID"  5  7 's/!integer/ integer/g'
     sed_range digda_grid.f90 "$ANCH_GRID" 10 12 's/ integer/!integer/g'
     sed_range digda_grid.f90 "$ANCH_GRID" 20 22 's/ integer/!integer/g'

    ## Mkinclude
    #sed -i -e "s/digda_256192/digda_128064/g" Mkinclude
    #sed -i -e "s/digda_21601080/digda_128064/g" Mkinclude
     ;;

   "256x192")
     # digda_mpi.f90
     sed_range digda_mpi.f90 "$ANCH_MPI" 0 13 's/ integer/!integer/g'
     sed_line  digda_mpi.f90 "$ANCH_MPI" 6    's/!integer/ integer/g'
     ncpu="48"

     # digda_grid.f90
     sed_range digda_grid.f90 "$ANCH_GRID"  5  7 's/ integer/!integer/g'
     sed_range digda_grid.f90 "$ANCH_GRID" 10 12 's/!integer/ integer/g'
     sed_range digda_grid.f90 "$ANCH_GRID" 20 22 's/ integer/!integer/g'

    #sed -i -e "s/digda_128064/digda_256192/g" Mkinclude
    #sed -i -e "s/digda_21601080/digda_256192/g" Mkinclude
     ;;

   "2160x1080")
     # digda_mpi.f90
     sed_range digda_mpi.f90 "$ANCH_MPI" 0 13 's/ integer/!integer/g'
     sed_line  digda_mpi.f90 "$ANCH_MPI" 8    's/!integer/ integer/g'
     ncpu="90"

     # digda_grid.f90
     sed_range digda_grid.f90 "$ANCH_GRID"  5  7 's/ integer/!integer/g'
     sed_range digda_grid.f90 "$ANCH_GRID" 10 12 's/ integer/!integer/g'
     sed_range digda_grid.f90 "$ANCH_GRID" 20 22 's/!integer/ integer/g'

    #sed -i -e "s/digda_128064/digda_21601080/g" Mkinclude
    #sed -i -e "s/digda_256192/digda_2160108/g" Mkinclude
     ;;

   *)
     echo "Unknown resolution: $resolution" >&2
     exit 1
     ;;
 esac

 # Update Mkinclude line with resolution and ncpu
 basemk=$(anchor_line Mkinclude "$ANCH_MK")
 sed -i -e "${basemk}s/digda/digda_${resolution}.${ncpu}${ciso}/g" Mkinclude

 echo "done: resolution=$resolution ncpu=$ncpu"

 ##############################

#make clean
#make > make.log 2>&1 &
#MAKE_PID=$!
#tail -f make.log &
#TAIL_PID=$!
#wait $MAKE_PID
#kill $TAIL_PID

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

 sed -i -e "${basemk}s/digda_${resolution}.${ncpu}${ciso}/digda/g" Mkinclude

 ##############################
