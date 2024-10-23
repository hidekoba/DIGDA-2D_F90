#!/bin/tcsh

 set dd = "/work/gi55/c24116/digda"
 set pd = "$dd/DIGDA-2D_F90/OUT"
 set ed = "$dd/DIGDA-2D_F90/OUT/1p4d/150000"
 set ed = "$dd/DIGDA-2D_F90/OUT/rccar_p4/1p4d/30000"
 set ed = "$dd/DIGDA-2D_F90/OUT/rccar_p4/1p4d/150000"

 ln -s $dd/ocnsed_coupling/cdo_gridfile/gridfile_2160x1080x1_lon_from0.txt gridfile_2160x1080x1.txt
 cdo remapbil,gridfile_2160x1080x1.txt $ed/SEDTRR.nc updated_SEDTRR.nc
 rm gridfile_2160x1080x1.txt
