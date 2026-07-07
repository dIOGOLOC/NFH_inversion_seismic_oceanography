# ===========================================
# Function to estimate parameters via TEOS-10
# ===========================================

import gsw
import xarray as xr
import numpy as np

def calc_pressure(depth, lat):
    
    """
    Calculate sea pressure from depth using TEOS-10.

    Parameters
    ----------
    depth : float
        Depth below the sea surface in meters (positive down).
    lat : float
        Latitude in decimal degrees, used to account for variations
        in gravity with latitude.

    Returns
    -------
    pressure : float
        Sea pressure in decibars (dbar), estimated from depth using
        the TEOS-10 equation of state (GSW function p_from_z).

    Notes
    -----
    This function converts depth to sea pressure assuming the reference
    dynamic height and sea surface geopotential are zero.
    """
     
    return gsw.p_from_z(-1*float(depth), float(lat), geo_strf_dyn_height=0, sea_surface_geopotential=0)


def estimate_tEOS10_parameters(ds):
    
    """
    Estimate key seawater thermodynamic parameters using TEOS-10.

    This function processes an xarray.Dataset containing oceanographic
    variables and computes derived quantities following TEOS-10 standards
    using the Gibbs Seawater (GSW) Oceanographic Toolbox.

    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset containing at least the following variables:
        - depth (m)
        - latitude (degrees north)
        - longitude (degrees east)
        - so : Practical Salinity (unitless)
        - thetao : Potential Temperature (°C)

    Returns
    -------
    ds : xarray.Dataset
        The input dataset with the following additional computed variables:

        - pressure (dbar)
            Sea pressure calculated from depth.
        - absolute_salinity (g/kg)
            Absolute Salinity computed from Practical Salinity.
        - conservative_temperature (°C)
            Conservative Temperature derived from potential temperature.
        - density (kg/m³)
            In-situ density calculated from Absolute Salinity, Conservative
            Temperature, and pressure.
        - sound_speed (m/s)
            Speed of sound in seawater derived from TEOS-10 relationships.

    Notes
    -----
    - All calculations follow TEOS-10 standards using the GSW-Python library.
    - Dask parallelization and xarray's vectorized ufuncs are used for efficiency.
    - Global dataset attributes are updated to record processing metadata,
      including TEOS-10 version and processing date.
    """

    # 1. Pressure (p)
    
    p = xr.apply_ufunc(
        calc_pressure,
        ds.depth, ds.latitude,
        input_core_dims=[[], []],
        output_core_dims=[[]],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[np.float64]
    )
    
    ds['pressure'] = p
    ds.pressure.attrs = {
        'long_name': 'Sea Pressure',
        'units': 'dbar',
        'standard_name': 'sea_water_pressure',
        'description': 'Sea pressure calculated from depth using TEOS-10',
        'positive': 'down'
    }
    
    # 2. Absolute Salinity (SA)

    SA = xr.apply_ufunc(
        gsw.SA_from_SP,
        ds.so, ds.pressure, ds.longitude, ds.latitude,
        input_core_dims=[[], [], [], []],
        output_core_dims=[[]],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[np.float64]
    )
    
    ds['absolute_salinity'] = SA
    ds.absolute_salinity.attrs = {
        'long_name': 'Absolute Salinity',
        'units': 'g/kg',
        'standard_name': 'sea_water_absolute_salinity',
        'description': 'Absolute Salinity calculated from Practical Salinity using TEOS-10'
    }
    
    # 3. Conservative Temperature (CT)
    
    CT = xr.apply_ufunc(
        gsw.CT_from_pt,
        ds.absolute_salinity, ds.thetao,
        input_core_dims=[[], []],
        output_core_dims=[[]],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[np.float64]
    )
    
    ds['conservative_temperature'] = CT
    ds.conservative_temperature.attrs = {
        'long_name': 'Conservative Temperature',
        'units': '°C',
        'standard_name': 'sea_water_conservative_temperature',
        'description': 'Conservative Temperature calculated from potential temperature using TEOS-10'
    }
    
    # 4. Density in-situ (rho)

    rho = xr.apply_ufunc(
        gsw.rho,
        ds.absolute_salinity, ds.conservative_temperature, ds.pressure,
        input_core_dims=[[], [], []],
        output_core_dims=[[]],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[np.float64]
    )
    
    ds['density'] = rho
    ds.density.attrs = {
        'long_name': 'In-situ Density',
        'units': 'kg/m³',
        'standard_name': 'sea_water_density',
        'description': 'In-situ density calculated using TEOS-10'
    }
    
    # 5. Sound Speed (ss)
    
    ss = xr.apply_ufunc(
        gsw.sound_speed,
        ds.absolute_salinity, ds.conservative_temperature, ds.pressure,
        input_core_dims=[[], [], []],
        output_core_dims=[[]],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[np.float64]
    )
    
    ds['sound_speed'] = ss
    ds.sound_speed.attrs = {
        'long_name': 'Sound Speed',
        'units': 'm/s',
        'standard_name': 'speed_of_sound_in_sea_water',
        'description': 'Sound speed in seawater calculated using TEOS-10'
    }
    
    
    # Adicionar atributos globais sobre o processamento
    ds.attrs['TEOS10_processing'] = 'Processed using GSW-Python and TEOS-10 standards'
    ds.attrs['processing_date'] = str(np.datetime64('now'))
    ds.attrs['GSW_version'] = gsw.__version__

    return ds