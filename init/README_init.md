# DIGDA initial-condition generator (standard / ciso)

Run from any directory (paths are resolved relative to the script):

```bash
bash init/mkinit_template_dummy.sh 128x64 standard
bash init/mkinit_template_dummy.sh 128x64 ciso
bash init/mkinit_template_dummy.sh 256x192 ciso
```

Select from `128x64`, `256x192`, and `2160x1080`. The second argument defaults to `standard`.

The `standard` mode generates the ten standard sediment variables; `ciso` adds
`C13-OrgM`, `C13-CaCO3`, `C13-DIC`, `C14-OrgM`, `C14-CaCO3`, and `C14-DIC`.
The ciso values reproduce the user-supplied example: zeros for particulate
isotope variables, 25 mmol m-3 for C13-DIC, and 2.6e-9 mmol m-3 for C14-DIC.
These numbers are **example initial concentrations**, not values automatically
calculated from isotope ratios or boundary forcing.

Requirements: Python 3, numpy, netCDF4, and NCO (`ncap2`, `ncks`, `ncwa`).
Five bottom-water files (`bnit.nc`, `bdic.nc`, `balk.nc`, `bsil.nc`, `bdo2.nc`)
are required under `Bound/data_digda<TAG>/`, or set `BOUND_DIR`.
Set `OUTPUT_DIR` to choose output location. Existing results are not overwritten.

`standard`: `SEDTRR0_<TAG>.nc`.
`ciso`: `SEDTRR0_<TAG>_ciso.nc` (rename/copy for a run script expecting another filename).

The five standard dissolved variables are initialized from boundary NetCDF
fields by NCO in either mode. Isotope variables use the listed uniform
example constants; they are **not** read from `bc13.nc` or `bc14.nc`.
Ensure the model executable and its `calc13`/`calc14` namelist switches support
the desired isotope configuration, and verify the restart variable names.

Runtime validation with actual boundary NetCDF files has not been performed.
