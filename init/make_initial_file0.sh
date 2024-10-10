#!/bin/bash

#dir="/work/01/gi55/c24116/digda/DIGDA-2D_F90/OUT/output_256192_500y_240722"
 dir="/work/01/gi55/c24116/digda/DIGDA-2D_F90/OUT/output_21601080_500y_240720"

 base_file="scla.nc"
 cp "$dir/$base_file" SEDTRR.nc

 file_list=("spom.nc" "scar.nc" "sopa.nc" "sage.nc" "cnit.nc" "cdic.nc" "calk.nc" "csil.nc" "cdo2.nc")
#file_list=("spom.nc" "scar.nc" "sopa.nc" "sage.nc" "cnit.nc" "cdic.nc" "calk.nc" "csil.nc" "cdo2.nc" "s13pom.nc" "s13car.nc" "c13dic.nc" "s14pom.nc" $dir/s14car.nc" "c14dic.nc")

 for file in "${file_list[@]}"; do
   ncks -A "$dir/$file" SEDTRR.nc
 done

#ncks -A $dir/scla.nc \
#        $dir/spom.nc \
#        $dir/scar.nc \
#        $dir/sopa.nc \
#        $dir/sage.nc \
#        $dir/cnit.nc \
#        $dir/cdic.nc \
#        $dir/calk.nc \
#        $dir/csil.nc \
#        $dir/cdo2.nc \
