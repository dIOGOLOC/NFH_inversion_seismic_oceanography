# Functions 

## modeling 

import numpy as np

def calculate_parameters(densities):
    '''
    Function to estimate the P and S velocities in function of initial density profile

    == P-wave and S-wave velocities to predict density == 
    
    Gardner, G.H.F., Gardner, L.W., and Gregory, A.R., 1974, Formation velocity and density – the
    diagnostic basics for stratigraphic traps: Geophysics, 39, 770-780.
    
    The purpose of this paper is to set forth certain relationships between rock physical properties,
    rock composition, and environmental conditions which have been established through extensive laboratory
    and field experimentation together with theoretical considerations. The literature on the subject is vast. 
    We are concerned primarily with seismic P-wave velocity and density of different types of sedimentary rocks in different
    environments. 
    
    Gardner et al. (1974) conducted a series of empirical studies and determined the following relationship between velocity and density:

    rho = a*Vp**(alpha)

    here \rho is in g/cm3, a is 0.31 when V is in m/s and alpha is 0.25.    
    

    == Vp/Vs ratio == 

    Vp/Vs velocity ratio is a strong function of:
    - water saturation, 
    - porosity, 
    - crack intensity, and
    - clay content. 
    
    Recently, it became important factor to study underground properties. 
    Vp/Vs velocity ratio has the variation interval as 1.45 to 8 and 
    it have been used as a lithological indicators in studies of soil amplification 
    and soil classification, acquifers and hydrocarbon reservoirs.

    For this work, we use the Vp/Vs equal square root of 3.


    Parameters:
    - densities (list) – Density profile in g/cm³.    
    
    Returns:
    vp, vs velocities in m/s (numpy.array)    
    '''
    
    vp = (densities/0.31)**4
    vs = vp / np.sqrt(3)
    
    return vp,vs


def calculate_parameters_from_vs(vel_s):
    '''
    Function to estimate the P velocities and Densities in function of initial Shear wave velocity
    
    == P-wave and S-wave velocities to predict density == 
    
    Gardner, G.H.F., Gardner, L.W., and Gregory, A.R., 1974, Formation velocity and density – the
    diagnostic basics for stratigraphic traps: Geophysics, 39, 770-780.
    
    The purpose of this paper is to set forth certain relationships between rock physical properties,
    rock composition, and environmental conditions which have been established through extensive laboratory
    and field experimentation together with theoretical considerations. The literature on the subject is vast. 
    We are concerned primarily with seismic P-wave velocity and density of different types of sedimentary rocks in different
    environments. 
    
    Gardner et al. (1974) conducted a series of empirical studies and determined the following relationship between velocity and density:

    rho = a*Vp**(alpha)

    here rho is in g/cm3, a is 0.31 when V is in m/s and alpha is 0.25.    
    

    == Vp/Vs ratio == 

    Vp/Vs velocity ratio is a strong function of:
    - water saturation, 
    - porosity, 
    - crack intensity, and
    - clay content. 
    
    Recently, it became important factor to study underground properties. 
    Vp/Vs velocity ratio has the variation interval as 1.45 to 8 and 
    it have been used as a lithological indicators in studies of soil amplification 
    and soil classification, acquifers and hydrocarbon reservoirs.

    For this work, we use the Vp/Vs equal square root of 3.


    Parameters:
    - vs (list or numpy.array) – Velocities in m/s.    
    
    Returns:
    vp velocities in m/s and density in g/cm³ (numpy.array)    
    ''' 
    
    vp = vel_s * np.sqrt(3)

    densities = 0.31*(vp**(1/4))

    return vp,densities*1000


# Function to generate a complete seismic model with random variation at each point
def create_seismic_model(depth_ranges,density_ranges,num_layers,xi,yj):
    '''
    xi: numpy.ndarray region [west-east points]) with a specific spacing (100 units in this case).
    yj: numpy.ndarray region [south-north points]) with a specific spacing (100 units in this case).
    num_layers: Number of layers of the model 
    '''
    
    # Initialize lists to store the variations of each parameter for each receiver
    x_l = []
    y_l = []
    depths = []
    vp = []
    vs = []
    densities = []
    formation = []
    
    for l in range(num_layers):
        x_l.append(xi)
        y_l.append(yj)
        
        # Generate random variations for each layer
        depths.append(np.random.uniform(*depth_ranges[l], size=xi.shape))
        densi = np.random.uniform(*density_ranges[l], size=yj.shape)
        densities.append(densi)
            
        # Calculate the parameters individually for each point
        vp.append(calculate_parameters(densi)[0])
        vs.append(calculate_parameters(densi)[1])

        formation.append('soil'+str(l+1))

    # Salvar modelo e dados
    model = {
        'xi':x_l,
        'yj':y_l,
        'vp': vp,
        'vs': vs,
        'depths': depths,
        'densities': densities,
        'num_layers': num_layers,
        'formation': formation,
    }

    return model
    
