#!/bin/bash

 dd="/work/gi55/c24116/digda"
 pd="$dd/DIGDA-2D_F90/OUT"
 ed="$pd/1p4d/150000"
 ed="$pd/rccar_p4/1p4d/30000"
 ed="$pd/rccar_p4/1p4d/150000"

 # Define input files
 fiocar="$ed/fiocar.nc"
 fiopom="$ed/fiopom.nc"
 fiopal="$ed/fioopa.nc"

 # Output text files
 output_orgm="global_OrgM_burial.txt"
 output_caco3="global_CaCO3_burial.txt"
 output_sio2="global_SiO2_burial.txt"

 # Final combined output text file
 combined_output="global_burial_fluxes.txt"

 # Conversion factor from mMol m-3 * cm s-1 to mol/yr
 # 1 mMol m-3 * cm s-1 * (1e-2 m/cm) * (1e-3 mol/mMol) * (3.1536e7 s/yr)
 conversion_factor=315.36

 # Conversion from mol to Tmol (1 Tmol = 1e12 mol)
 mol_to_tmol=1e-12
 # Conversion from mol to PgC (1 mol C = 12 g C, 1 PgC = 1e15 g C)
 mol_to_pgc=12e-15

 # Function to calculate the global burial flux and convert units
 calculate_burial_flux() {
   input_file=$1
   variable_name=$2
   output_file=$3
   conversion_factor_pgc=$4

   # Calculate the grid cell area in m^2
   cdo gridarea "$input_file" grid_area.nc

   # Multiply burial flux by grid area
   cdo mul "$input_file" grid_area.nc tmp_burial_area.nc

   # Perform the sum over space (latitude and longitude)
   cdo fldsum tmp_burial_area.nc spatial_sum.nc

   # Multiply by conversion factor to convert units to mol/yr
   cdo mulc,$conversion_factor spatial_sum.nc converted_sum.nc

#### Extract the result and write to a text file
###cdo outputtab,value converted_sum.nc > "$output_file"

   # Extract the Tmol/yr result and write to a text file
   cdo mulc,$mol_to_tmol converted_sum.nc tmp_tmol.nc
   cdo outputtab,value tmp_tmol.nc > "$output_file"

   # If PgC conversion is requested
   if [ "$conversion_factor_pgc" = true ]; then
     cdo mulc,$mol_to_pgc converted_sum.nc tmp_pgc.nc
     cdo outputtab,value tmp_pgc.nc >> "$output_file"
   fi

   # Clean up temporary files
   rm -f grid_area.nc tmp_burial_area.nc spatial_sum.nc converted_sum.nc tmp_tmol.nc tmp_pgc.nc
 }

 # Calculate burial flux for each variable
 calculate_burial_flux "$fiopom" "OrgM_burial"  "$output_orgm"  true
 calculate_burial_flux "$fiocar" "CaCO3_burial" "$output_caco3" true
 calculate_burial_flux "$fiopal" "SiO2_burial"  "$output_sio2"  false

 # Combine the results into a single text file
 echo "OrgM_burial (Tmol/yr, PgC/yr):" > "$combined_output"
 cat "$output_orgm" >> "$combined_output"
 echo "" >> "$combined_output"

 echo "CaCO3_burial (Tmol/yr, PgC/yr):" >> "$combined_output"
 cat "$output_caco3" >> "$combined_output"
 echo "" >> "$combined_output"

 echo "SiO2_burial (Tmol/yr, PgC/yr):" >> "$combined_output"
 cat "$output_sio2" >> "$combined_output"

 # Clean up individual output text files
 rm -f "$output_caco3" "$output_orgm" "$output_sio2"

 echo "Global burial fluxes have been calculated and saved to text files."
