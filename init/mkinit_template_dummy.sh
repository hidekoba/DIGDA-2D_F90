#!/usr/bin/env bash
# Example initialization of DIGDA restart data.
# Run from init/: bash mkinit_template_dummy.sh 128x64
set -euo pipefail

resolution="${1:-256x192}"
mode="${2:-standard}"
case "$mode" in
  standard|ciso) ;;
  *) echo "Usage: $0 {128x64|256x192|2160x1080} [standard|ciso]" >&2; exit 2 ;;
esac
case "$resolution" in
  128x64)    tag="128064";   bound_name="data_digda128064" ;;
  256x192)   tag="256192";   bound_name="data_digda256192" ;;
  2160x1080) tag="21601080"; bound_name="data_digda21601080" ;;
  *) echo "Usage: $0 {128x64|256x192|2160x1080}" >&2; exit 2 ;;
esac

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bound_dir="${BOUND_DIR:-$here/../Bound/$bound_name}"
output_dir="${OUTPUT_DIR:-$here}"
mkdir -p "$output_dir"

for command in python3 ncap2 ncks ncwa; do
  command -v "$command" >/dev/null || { echo "Missing command: $command" >&2; exit 1; }
done
python3 -c 'import netCDF4, numpy' || { echo 'Install Python netCDF4 and numpy' >&2; exit 1; }
for name in bnit bdic balk bsil bdo2; do
  [[ -f "$bound_dir/$name.nc" ]] || { echo "Missing: $bound_dir/$name.nc" >&2; exit 1; }
done
[[ -f "$here/mkinit_template_${tag}.py" ]] || { echo "Missing Python template for $resolution" >&2; exit 1; }

# Work in a temporary directory so source files and previous outputs are protected.
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/digda_init.XXXXXXXX")"
trap 'rm -rf -- "$work_dir"' EXIT
cp "$here/mkinit_template_${tag}.py" "$work_dir/"
cp "$here/"{lon,lat,dep}_values_${tag}10.txt "$work_dir/"
cd "$work_dir"
if [[ "$mode" == "ciso" ]]; then
  python3 "mkinit_template_${tag}.py" --ciso
else
  python3 "mkinit_template_${tag}.py"
fi
cp template.nc temp.nc

# Expand each (time,lat,lon) bottom-water field into ten sediment layers.
# The source variable names and units must match the model's boundary data.
for spec in 'bnit:NIT:Phosphate' 'bdic:DIC:DIC' 'balk:ALK:Alkalinity' 'bsil:SIL:Silicate' 'bdo2:DO2:Oxygen'; do
  IFS=: read -r fname inputvar outputvar <<< "$spec"
  ncap2 -O -s "defdim(\"depth\",10); ${inputvar}_3D[\$time,\$depth,\$lat,\$lon]=${inputvar}(\$time,\$lat,\$lon)" "$bound_dir/$fname.nc" file3D.nc
  ncap2 -O -s 'time=int(time)' file3D.nc file3D_int.nc
  ncwa -O -a time file3D_int.nc file3D_int_zero_dim.nc
  ncks -A -v "${inputvar}_3D" file3D_int_zero_dim.nc temp.nc
  ncap2 -O -s "${outputvar}=${inputvar}_3D" temp.nc updated_template.nc
  mv updated_template.nc temp.nc
  rm -f file3D.nc file3D_int.nc file3D_int_zero_dim.nc
done
ncks -C -O -x -v NIT_3D,DIC_3D,ALK_3D,SIL_3D,DO2_3D temp.nc final.nc
# Never silently replace an earlier initialization file.
if [[ "$mode" == "ciso" ]]; then
  out="$output_dir/SEDTRR0_${tag}_ciso.nc"
else
  out="$output_dir/SEDTRR0_${tag}.nc"
fi
if [[ -e "$out" ]]; then
  echo "Output already exists: $out (move or remove it explicitly)" >&2
  exit 1
fi
mv final.nc "$out"
echo "Created: $out"
