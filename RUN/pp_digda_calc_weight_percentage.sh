 #!/bin/bash

 # Check if the directory is passed as an argument
 if [ "$#" -ne 1 ]; then
     echo "Usage: $0 directory"
     exit 1
 fi

 dir=$1

 # Check if the specified directory exists
 if [ ! -d "$dir" ]; then
     echo "Directory $dir does not exist."
     exit 1
 fi

 # Load the conda initialization script
 source /home/c24116/anaconda3/etc/profile.d/conda.sh

 # Activate an Anaconda environment
 conda activate /home/c24116/anaconda3/

 # Extract relevant variables from each file and combine into a single file for calculation
 rm combined.nc
 ###ncks -A -v Clay  $dir/scla.nc combined.nc
 ###ncks -A -v OrgM  $dir/spom.nc combined.nc
 ###ncks -A -v CaCO3 $dir/scar.nc combined.nc
 ###ncks -A -v SiO2  $dir/sopa.nc combined.nc
 ncks -A -v Clay  $dir/SEDTRR.nc combined.nc
 ncks -A -v OrgM  $dir/SEDTRR.nc combined.nc
 ncks -A -v CaCO3 $dir/SEDTRR.nc combined.nc
 ncks -A -v SiO2  $dir/SEDTRR.nc combined.nc


 # Calculate WT and handle cases where WT == 0
 ncap2 -O -s 'MOLPOM=12.0; MOLCAR=100.0; MOLOPA=60.0; WT=Clay+(OrgM*MOLPOM*1.0e-3)+(CaCO3*MOLCAR*1.0e-3)+(SiO2*MOLOPA*1.0e-3)' combined.nc wt.nc
 ncap2 -O -s 'WT=(WT==0)?1.0:WT' wt.nc wt.nc  # Avoid division by zero by setting WT to 1 where it was zero

 # Calculate PCLA
 ncap2 -O -s 'PCLA=(WT==1.0)?0.0:(Clay/WT)*1.0e2' wt.nc pcla.nc

 # Calculate PPOM
 ncap2 -O -s 'PPOM=(WT==1.0)?0.0:(OrgM*MOLPOM*1.0e-3/WT)*1.0e2' wt.nc ppom.nc

 # Calculate PCAR
 ncap2 -O -s 'PCAR=(WT==1.0)?0.0:(CaCO3*MOLCAR*1.0e-3/WT)*1.0e2' wt.nc pcar.nc

 # Calculate POPA
 ncap2 -O -s 'POPA=(WT==1.0)?0.0:(SiO2*MOLOPA*1.0e-3/WT)*1.0e2' wt.nc popa.nc

 # Delete unnecessary variables
 ncks -O -x -v WT,MOLCAR,MOLOPA,MOLPOM,Clay,OrgM,CaCO3,SiO2 pcla.nc pcla.nc
 ncks -O -x -v WT,MOLCAR,MOLOPA,MOLPOM,Clay,OrgM,CaCO3,SiO2 ppom.nc ppom.nc
 ncks -O -x -v WT,MOLCAR,MOLOPA,MOLPOM,Clay,OrgM,CaCO3,SiO2 pcar.nc pcar.nc
 ncks -O -x -v WT,MOLCAR,MOLOPA,MOLPOM,Clay,OrgM,CaCO3,SiO2 popa.nc popa.nc

 # Setting variable attributes
 ncatted -O -a _FillValue,PCLA,o,d,999. pcla.nc
 ncatted -O -a cell_methods,PCLA,o,c,"time: mean" pcla.nc
 ncatted -O -a missing_value,PCLA,o,d,999. pcla.nc
 ncatted -O -a original_name,PCLA,o,c,"Clay" pcla.nc
 ncatted -O -a units,PCLA,o,c,"ND" pcla.nc

 ncatted -O -a _FillValue,depth,o,d,999. pcla.nc
 ncatted -O -a standard_name,depth,o,c,"depth" pcla.nc
 ncatted -O -a long_name,depth,o,c,"depth" pcla.nc
 ncatted -O -a units,depth,o,c,"m" pcla.nc
 ncatted -O -a axis,depth,o,c,"Z" pcla.nc

 ncatted -O -a standard_name,lat,o,c,"latitude" pcla.nc
 ncatted -O -a long_name,lat,o,c,"latitude" pcla.nc
 ncatted -O -a units,lat,o,c,"degrees_north" pcla.nc
 ncatted -O -a axis,lat,o,c,"Y" pcla.nc

 ncatted -O -a standard_name,lon,o,c,"longitude" pcla.nc
 ncatted -O -a long_name,lon,o,c,"longitude" pcla.nc
 ncatted -O -a units,lon,o,c,"degrees_east" pcla.nc
 ncatted -O -a axis,lon,o,c,"X" pcla.nc

 # Setting variable attributes
 ncatted -O -a _FillValue,PPOM,o,d,999. ppom.nc
 ncatted -O -a cell_methods,PPOM,o,c,"time: mean" ppom.nc
 ncatted -O -a missing_value,PPOM,o,d,999. ppom.nc
 ncatted -O -a original_name,PPOM,o,c,"OrgM" ppom.nc
 ncatted -O -a units,PPOM,o,c,"ND" ppom.nc

 ncatted -O -a _FillValue,depth,o,d,999. ppom.nc
 ncatted -O -a standard_name,depth,o,c,"depth" ppom.nc
 ncatted -O -a long_name,depth,o,c,"depth" ppom.nc
 ncatted -O -a units,depth,o,c,"m" ppom.nc
 ncatted -O -a axis,depth,o,c,"Z" ppom.nc

 ncatted -O -a standard_name,lat,o,c,"latitude" ppom.nc
 ncatted -O -a long_name,lat,o,c,"latitude" ppom.nc
 ncatted -O -a units,lat,o,c,"degrees_north" ppom.nc
 ncatted -O -a axis,lat,o,c,"Y" ppom.nc

 ncatted -O -a standard_name,lon,o,c,"longitude" ppom.nc
 ncatted -O -a long_name,lon,o,c,"longitude" ppom.nc
 ncatted -O -a units,lon,o,c,"degrees_east" ppom.nc
 ncatted -O -a axis,lon,o,c,"X" ppom.nc

 # Setting variable attributes
 ncatted -O -a _FillValue,PCAR,o,d,999. pcar.nc
 ncatted -O -a cell_methods,PCAR,o,c,"time: mean" pcar.nc
 ncatted -O -a missing_value,PCAR,o,d,999. pcar.nc
 ncatted -O -a original_name,PCAR,o,c,"CaCO3" pcar.nc
 ncatted -O -a units,PCAR,o,c,"ND" pcar.nc

 ncatted -O -a _FillValue,depth,o,d,999. pcar.nc
 ncatted -O -a standard_name,depth,o,c,"depth" pcar.nc
 ncatted -O -a long_name,depth,o,c,"depth" pcar.nc
 ncatted -O -a units,depth,o,c,"m" pcar.nc
 ncatted -O -a axis,depth,o,c,"Z" pcar.nc

 ncatted -O -a standard_name,lat,o,c,"latitude" pcar.nc
 ncatted -O -a long_name,lat,o,c,"latitude" pcar.nc
 ncatted -O -a units,lat,o,c,"degrees_north" pcar.nc
 ncatted -O -a axis,lat,o,c,"Y" pcar.nc

 ncatted -O -a standard_name,lon,o,c,"longitude" pcar.nc
 ncatted -O -a long_name,lon,o,c,"longitude" pcar.nc
 ncatted -O -a units,lon,o,c,"degrees_east" pcar.nc
 ncatted -O -a axis,lon,o,c,"X" pcar.nc

 # Setting variable attributes
 ncatted -O -a _FillValue,POPA,o,d,999. popa.nc
 ncatted -O -a cell_methods,POPA,o,c,"time: mean" popa.nc
 ncatted -O -a missing_value,POPA,o,d,999. popa.nc
 ncatted -O -a original_name,POPA,o,c,"SiO2" popa.nc
 ncatted -O -a units,POPA,o,c,"ND" popa.nc

 ncatted -O -a _FillValue,depth,o,d,999. popa.nc
 ncatted -O -a standard_name,depth,o,c,"depth" popa.nc
 ncatted -O -a long_name,depth,o,c,"depth" popa.nc
 ncatted -O -a units,depth,o,c,"m" popa.nc
 ncatted -O -a axis,depth,o,c,"Z" popa.nc

 ncatted -O -a standard_name,lat,o,c,"latitude" popa.nc
 ncatted -O -a long_name,lat,o,c,"latitude" popa.nc
 ncatted -O -a units,lat,o,c,"degrees_north" popa.nc
 ncatted -O -a axis,lat,o,c,"Y" popa.nc

 ncatted -O -a standard_name,lon,o,c,"longitude" popa.nc
 ncatted -O -a long_name,lon,o,c,"longitude" popa.nc
 ncatted -O -a units,lon,o,c,"degrees_east" popa.nc
 ncatted -O -a axis,lon,o,c,"X" popa.nc

 # Calculate global average for each variable and output the result to the terminal
 output_file="wtp_global_averages.txt"
 ncap2 -O -s 'PCLA_avg=avg(PCLA)' pcla.nc -o tmp.nc ; PCLA_avg=$(ncks -v PCLA_avg tmp.nc | grep "PCLA_avg =" | awk '{print $3}')
 ncap2 -O -s 'PPOM_avg=avg(PPOM)' ppom.nc -o tmp.nc ; PPOM_avg=$(ncks -v PPOM_avg tmp.nc | grep "PPOM_avg =" | awk '{print $3}')
 ncap2 -O -s 'PCAR_avg=avg(PCAR)' pcar.nc -o tmp.nc ; PCAR_avg=$(ncks -v PCAR_avg tmp.nc | grep "PCAR_avg =" | awk '{print $3}')
 ncap2 -O -s 'POPA_avg=avg(POPA)' popa.nc -o tmp.nc ; POPA_avg=$(ncks -v POPA_avg tmp.nc | grep "POPA_avg =" | awk '{print $3}')

 # Output the global average values
 echo "Global Average of PCLA: $PCLA_avg" >  $output_file
 echo "Global Average of PPOM: $PPOM_avg" >> $output_file
 echo "Global Average of PCAR: $PCAR_avg" >> $output_file
 echo "Global Average of POPA: $POPA_avg" >> $output_file

 # Move the generated files to that directory
 rm combined.nc wt.nc
 mv pcla.nc ppom.nc pcar.nc popa.nc $dir
 mv $output_file $dir

 # Deactivate the conda environment after execution
 conda deactivate
