#!/bin/bash

 # Load the conda initialization script
 source /home/c24116/anaconda3/etc/profile.d/conda.sh
 # Activate an Anaconda environment
 conda activate /home/c24116/anaconda3/

 ##################################################
 # make initial restart file
 #   256x192
 ##################################################

 python mkinit_template_256192.py

 cp template.nc temp.nc

 ncap2 -O -s 'defdim("depth",10); NIT_3D[$time,$depth,$lat,$lon]=NIT($time,$lat,$lon)' ../Bound/data_digda256192/bnit.nc file3D.nc
 ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ncks -A -v NIT_3D file3D_int_zero_dim.nc temp.nc
 ncap2 -O -s 'Phosphate=NIT_3D'  temp.nc updated_template.nc
 rm file3D.nc file3D_int.nc file3D_int_zero_dim.nc

 ncap2 -O -s 'defdim("depth",10); DIC_3D[$time,$depth,$lat,$lon]=DIC($time,$lat,$lon)' ../Bound/data_digda256192/bdic.nc file3D.nc
 ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ncks -A -v DIC_3D file3D_int_zero_dim.nc temp.nc
 ncap2 -O -s 'DIC=DIC_3D'  temp.nc updated_template.nc

 ncap2 -O -s 'defdim("depth",10); ALK_3D[$time,$depth,$lat,$lon]=ALK($time,$lat,$lon)' ../Bound/data_digda256192/balk.nc file3D.nc
 ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ncks -A -v ALK_3D file3D_int_zero_dim.nc temp.nc
 ncap2 -O -s 'Alkalinity=ALK_3D'  temp.nc updated_template.nc

 ncap2 -O -s 'defdim("depth",10); SIL_3D[$time,$depth,$lat,$lon]=SIL($time,$lat,$lon)' ../Bound/data_digda256192/bsil.nc file3D.nc
 ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ncks -A -v SIL_3D file3D_int_zero_dim.nc temp.nc
 ncap2 -O -s 'Silicate=SIL_3D'  temp.nc updated_template.nc

 ncap2 -O -s 'defdim("depth",10); DO2_3D[$time,$depth,$lat,$lon]=DO2($time,$lat,$lon)' ../Bound/data_digda256192/bdo2.nc file3D.nc
 ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ncks -A -v DO2_3D file3D_int_zero_dim.nc temp.nc
 ncap2 -O -s 'Oxygen=DO2_3D'  temp.nc updated_template.nc

 ###ncap2 -O -s 'defdim("depth",10); C13_3D[$time,$depth,$lat,$lon]=C13($time,$lat,$lon)' ../Bound/data_digda256192/bc13.nc file3D.nc
 ###ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ###ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ###ncks -A -v C13_3D file3D_int_zero_dim.nc temp.nc
 ###ncap2 -O -s '13CDIC=C13_3D'  temp.nc updated_template.nc

 ###ncap2 -O -s 'defdim("depth",10); C14_3D[$time,$depth,$lat,$lon]=C14($time,$lat,$lon)' ../Bound/data_digda256192/bc14.nc file3D.nc
 ###ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ###ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ###ncks -A -v C14_3D file3D_int_zero_dim.nc temp.nc
 ###ncap2 -O -s '14CDIC=C14_3D'  temp.nc updated_template.nc

 rm temp.nc ; mv updated_template.nc template.nc
 ncks -C -O -x -v NIT_3D,DIC_3D,ALK_3D,SIL_3D,DO2_3D template.nc output.nc ; mv output.nc template.nc

 mv template.nc SEDTRR0_256192.nc

###################################################
## make initial restart file
##   2160x1080
###################################################

 python mkinit_template_21601080.py

 cp template.nc temp.nc

 ncap2 -O -s 'defdim("depth",10); NIT_3D[$time,$depth,$lat,$lon]=NIT($time,$lat,$lon)' ../Bound/data_digda21601080/bnit.nc file3D.nc
 ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ncks -A -v NIT_3D file3D_int_zero_dim.nc temp.nc
 ncap2 -O -s 'Phosphate=NIT_3D'  temp.nc updated_template.nc
 rm file3D.nc file3D_int.nc file3D_int_zero_dim.nc

 ncap2 -O -s 'defdim("depth",10); DIC_3D[$time,$depth,$lat,$lon]=DIC($time,$lat,$lon)' ../Bound/data_digda21601080/bdic.nc file3D.nc
 ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ncks -A -v DIC_3D file3D_int_zero_dim.nc temp.nc
 ncap2 -O -s 'DIC=DIC_3D'  temp.nc updated_template.nc

 ncap2 -O -s 'defdim("depth",10); ALK_3D[$time,$depth,$lat,$lon]=ALK($time,$lat,$lon)' ../Bound/data_digda21601080/balk.nc file3D.nc
 ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ncks -A -v ALK_3D file3D_int_zero_dim.nc temp.nc
 ncap2 -O -s 'Alkalinity=ALK_3D'  temp.nc updated_template.nc

 ncap2 -O -s 'defdim("depth",10); SIL_3D[$time,$depth,$lat,$lon]=SIL($time,$lat,$lon)' ../Bound/data_digda21601080/bsil.nc file3D.nc
 ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ncks -A -v SIL_3D file3D_int_zero_dim.nc temp.nc
 ncap2 -O -s 'Silicate=SIL_3D'  temp.nc updated_template.nc

 ncap2 -O -s 'defdim("depth",10); DO2_3D[$time,$depth,$lat,$lon]=DO2($time,$lat,$lon)' ../Bound/data_digda21601080/bdo2.nc file3D.nc
 ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ncks -A -v DO2_3D file3D_int_zero_dim.nc temp.nc
 ncap2 -O -s 'Oxygen=DO2_3D'  temp.nc updated_template.nc

 ###ncap2 -O -s 'defdim("depth",10); C13_3D[$time,$depth,$lat,$lon]=C13($time,$lat,$lon)' ../Bound/data_digda21601080/bc13.nc file3D.nc
 ###ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ###ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ###ncks -A -v C13_3D file3D_int_zero_dim.nc temp.nc
 ###ncap2 -O -s '13CDIC=C13_3D'  temp.nc updated_template.nc

 ###ncap2 -O -s 'defdim("depth",10); C14_3D[$time,$depth,$lat,$lon]=C14($time,$lat,$lon)' ../Bound/data_digda21601080/bc14.nc file3D.nc
 ###ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
 ###ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
 ###ncks -A -v C14_3D file3D_int_zero_dim.nc temp.nc
 ###ncap2 -O -s '14CDIC=C14_3D'  temp.nc updated_template.nc

 rm temp.nc ; mv updated_template.nc template.nc
 ncks -C -O -x -v NIT_3D,DIC_3D,ALK_3D,SIL_3D,DO2_3D template.nc output.nc ; mv output.nc template.nc

 mv template.nc SEDTRR0_21601080.nc

 # Deactivate the conda environment after execution
 conda deactivate
