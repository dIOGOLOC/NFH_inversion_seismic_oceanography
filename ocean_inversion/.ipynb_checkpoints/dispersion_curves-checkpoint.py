# Functions 

## Dispersion curve estimative 

import numpy as np
from disba import PhaseDispersion,EigenFunction,PhaseSensitivity
import pandas as pd

from CODES.modeling import create_seismic_model,calculate_parameters,calculate_parameters_from_vs

def create_velocity_model_from_profile(model_profile):
    '''
    Parameters
    ----------
    Density profile model: 1-D numpy.array
        [thickness,velocity_p,velocity_s,density]
            - Layer thickness (in km).
            - Layer P-wave velocity (in km/s).
            - Layer S-wave velocity (in km/s).
            - Layer density (in g/cm3).     
    '''
  

    dens_values_unique,dens_index_unique = np.unique(model_profile,return_index=True)

    dens_index_unique = sorted(dens_index_unique)

    depth_ = np.linspace(0,-2.,len(model_profile)+1)

    thickness_ = [depth_[i] for i in dens_index_unique]
    thickness_.append(-2.0)

    thickness = np.diff(thickness_)*(-1)

    dens = model_profile[dens_index_unique]

    vp, vs = calculate_parameters(dens)

    vmodel = []
    for idx in range(len(dens)):

        # thickness, Vp, Vs, density
        # km, km/s, km/s, g/cm3        

        if not idx == len(dens) - 1:
            vmodel.append([thickness[idx]/1000,vp[idx]/1000,vs[idx]/1000,dens[idx]/1000])
        else: 
            vmodel.append([0,vp[idx]/1000,vs[idx]/1000,dens[idx]/1000])
   
    velocity_model = np.array(vmodel)    

    return velocity_model

# -----------------------------------------------------------

def create_velocity_model_from_profile_vs(model_profile):
    '''
    Parameters
    ----------
    Velocity model: 1-D numpy.array
        [thickness,velocity_p,velocity_s,density]
            - Layer thickness (in km).
            - Layer P-wave velocity (in km/s).
            - Layer S-wave velocity (in km/s).
            - Layer density (in g/cm3).     
    '''
    
    vmodel = []
    for i, (thickness, vs) in enumerate(zip(*model_profile)):
 
        vp, dens = calculate_parameters_from_vs(vs)
        
        # thickness, Vp, Vs, density
        # km, km/s, km/s, g/cm3
        if not i == len(model_profile[0]) - 1:
            vmodel.append([thickness/1000,vp/1000,vs/1000,dens/1000])
        else: 
            vmodel.append([00,vp/1000,vs/1000,dens/1000])


   
    velocity_model = np.array(vmodel)    

    return velocity_model
    
# -----------------------------------------------------------

def estimate_disp_from_velocity_model(vel_mol,number_samples=100,algorithm_str='dunkin'):
    '''
    Calculate phase velocities for input velocity model.

    Parameters
    ----------
    Velocity model: 2-D numpy.array
        [thickness,velocity_p,velocity_s,density]
            - Layer thickness (in km).
            - Layer P-wave velocity (in km/s).
            - Layer S-wave velocity (in km/s).
            - Layer density (in g/cm3).
            
    ----
    Example of velocity model by DISBA (https://keurfonluu.github.io/disba/)
    
    # Velocity model
    # thickness, Vp, Vs, density
    # km, km/s, km/s, g/cm3
    velocity_model = np.array([
	   [10.0, 7.00, 3.50, 2.00],
	   [10.0, 6.80, 3.40, 2.00],
	   [10.0, 7.00, 3.50, 2.00],
	   [10.0, 7.60, 3.80, 2.00],
	   [10.0, 8.40, 4.20, 2.00],
	   [10.0, 9.00, 4.50, 2.00],
	   [10.0, 9.40, 4.70, 2.00],
	   [10.0, 9.60, 4.80, 2.00],
	   [10.0, 9.50, 4.75, 2.00],])

	# Periods must be sorted starting with low periods
	t = np.logspace(0.0, 3.0, 100)

	# Compute the 3 first Rayleigh- and Love- wave modal dispersion curves
	# Fundamental mode corresponds to mode 0
	pd = PhaseDispersion(*velocity_model.T)
	cpr = [pd(t, mode=i, wave="rayleigh") for i in range(3)]
	cpl = [pd(t, mode=i, wave="love") for i in range(3)]
	'''

    # Periods must be sorted starting with low periods
    hz = np.linspace(10, 100.0, number_samples) # Hertz
    
    t = 1/hz[::-1] # Hertz to seconds
    
    # Fundamental mode corresponds to mode 0
    pdisp = PhaseDispersion(*vel_mol.T,dc=0.0001,algorithm=algorithm_str)
    
    # Fundamental mode corresponds to mode 0
    cpr = pdisp(t, mode=0, wave="rayleigh")

    return cpr    

# -----------------------------------------------------------

def compute_dispersion(row,vs_col='mean_vs',depth_col='mean_depth'):
    
    # Extract Vs and thickness from the row
    Vs = row[vs_col]
    thick = row[depth_col]
    
    # Create velocity model and estimate dispersion curve
    simulated_velocity_model = create_velocity_model_from_profile_vs([thick, Vs])
    simulated_cpr = estimate_disp_from_velocity_model(simulated_velocity_model)
    
    # Extract simulated values
    simulated_dispersion = simulated_cpr.velocity * 1000  # Convert km/s to m/s
    simulated_frequency = 1/simulated_cpr.period         # Already in seconds (assuming)
    
    return pd.Series({
        'simulated_dispersion': simulated_dispersion,
        'simulated_frequency': simulated_frequency
    })