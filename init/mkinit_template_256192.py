import argparse
import netCDF4 as nc
import numpy as np
from datetime import datetime

parser = argparse.ArgumentParser(description="Create a DIGDA initial NetCDF template")
parser.add_argument("--ciso", action="store_true", help="Include carbon isotope tracers")
args = parser.parse_args()

# Path for the output NetCDF file
output_file = 'template.nc'

lon_values_file = 'lon_values_25619210.txt'
lat_values_file = 'lat_values_25619210.txt'
dep_values_file = 'dep_values_25619210.txt'

# Function to read lat values from a file
def read_dim_values(file_path):
    with open(file_path, 'r') as f:
        dim_values = np.array([float(x) for x in f.read().split(',')])
    return dim_values

# Create a new NetCDF file
with nc.Dataset(output_file, 'w', format='NETCDF4') as ds:
    # Define dimensions (lon, lat, depth, time)
    lon_size = 256
    lat_size = 192
    depth_size = 10
    time_size = None  # Unlimited dimension for time

    # Create dimensions in the NetCDF file
    ds.createDimension('lon', lon_size)
    ds.createDimension('lat', lat_size)
    ds.createDimension('depth', depth_size)
    ds.createDimension('time', time_size)  # Unlimited dimension

    # Create variables for dimensions
    lon = ds.createVariable('lon', 'f8', ('lon',))
    lat = ds.createVariable('lat', 'f8', ('lat',))
    depth = ds.createVariable('depth', 'f8', ('depth',))
    time = ds.createVariable('time', 'i4', ('time',))

    # Create a data variable (e.g., scla)
    scla = ds.createVariable('Clay','f8', ('time', 'depth', 'lat', 'lon'), fill_value=-9999.0)
    spom = ds.createVariable('OrgM', 'f8', ('time', 'depth', 'lat', 'lon'), fill_value=-9999.0)
    scar = ds.createVariable('CaCO3', 'f8', ('time', 'depth', 'lat', 'lon'), fill_value=-9999.0)
    sopa = ds.createVariable('SiO2', 'f8', ('time', 'depth', 'lat', 'lon'), fill_value=-9999.0)
    sage = ds.createVariable('Agetracer', 'f8', ('time', 'depth', 'lat', 'lon'), fill_value=-9999.0)
    cnit = ds.createVariable('Phosphate', 'f8', ('time', 'depth', 'lat', 'lon'), fill_value=-9999.0)
    cdic = ds.createVariable('DIC', 'f8', ('time', 'depth', 'lat', 'lon'), fill_value=-9999.0)
    calk = ds.createVariable('Alkalinity', 'f8', ('time', 'depth', 'lat', 'lon'), fill_value=-9999.0)
    csil = ds.createVariable('Silicate', 'f8', ('time', 'depth', 'lat', 'lon'), fill_value=-9999.0)
    cdo2 = ds.createVariable('Oxygen', 'f8', ('time', 'depth', 'lat', 'lon'), fill_value=-9999.0)

    if args.ciso:
        isotope_names = ('C13-OrgM', 'C13-CaCO3', 'C13-DIC',
                         'C14-OrgM', 'C14-CaCO3', 'C14-DIC')
        isotopes = {}
        for name in isotope_names:
            variable = ds.createVariable(name, 'f8',
                                         ('time', 'depth', 'lat', 'lon'),
                                         fill_value=-9999.0)
            variable.units = 'mmol m-3'
            variable.long_name = name
            isotopes[name] = variable

    # Set attributes for the variables
    time.units = 'year'
    time.calendar = 'proleptic_gregorian'
    scla.units = 'g m-3'
    spom.units = 'mmol m-3'
    scar.units = 'mmol m-3'
    sopa.units = 'mmol m-3'
    sage.units = 'mmol m-3'
    cnit.units = 'mmol m-3'
    cdic.units = 'mmol m-3'
    calk.units = 'mmol m-3'
    csil.units = 'mmol m-3'
    cdo2.units = 'mmol m-3'
    scla.long_name = 'Clay'
    spom.long_name = 'OrgM'
    scar.long_name = 'CaCO3'
    sopa.long_name = 'SiO2'
    sage.long_name = 'Agetracer'
    cnit.long_name = 'Phosphate'
    cdic.long_name = 'DIC'
    calk.long_name = 'Alkalinity'
    csil.long_name = 'Silicate'
    cdo2.long_name = 'Oxygen'

    # Define values for latitude, longitude, and depth
    lon_values = read_dim_values(lon_values_file)
    lat_values = read_dim_values(lat_values_file)
    depth_values = read_dim_values(dep_values_file)

    # Assign the values to the corresponding variables
    lon[:] = lon_values
    lat[:] = lat_values
    depth[:] = depth_values

    # Set time data (for example, using 0 as the first timestep)
    time[0] = 0

    # Initialize the scla variable with a default value
    scla[0, :, :, :] = 2600000.0
    spom[0, :, :, :] = 0.0
    scar[0, :, :, :] = 0.0
    sopa[0, :, :, :] = 0.0
    sage[0, :, :, :] = 0.0
    cnit[0, :, :, :] = 0.0
    cdic[0, :, :, :] = 0.0
    calk[0, :, :, :] = 0.0
    csil[0, :, :, :] = 0.0
    cdo2[0, :, :, :] = 0.0

    if args.ciso:
        # Example carbon isotope initial concentrations, as supplied in
        # the ciso initialization template. Adjust for your experiment.
        initial_isotopes = {
            'C13-OrgM': 0.0,
            'C13-CaCO3': 0.0,
            'C13-DIC': 25.0,
            'C14-OrgM': 0.0,
            'C14-CaCO3': 0.0,
            'C14-DIC': 2.6e-9,
        }
        for name, value in initial_isotopes.items():
            isotopes[name][0, :, :, :] = value

    # Set global attributes for the NetCDF file
    ds.description = 'NetCDF example created with Python'
    ds.history = f'Created on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

# Print success message
print(f'The NetCDF file {output_file} has been created successfully.')
