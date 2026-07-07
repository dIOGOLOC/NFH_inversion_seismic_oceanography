# Functions 

## Dispersion curve estimative 

import numpy as np
import pandas as pd
from itertools import accumulate
import multiprocessing


from deap import base, creator, tools, algorithms

from CODES.modeling import calculate_parameters_from_vs,create_velocity_model_from_profile_vs
from CODES.dispersion_curves import estimate_disp_from_velocity_model

from parameters_py.config import (
<<<<<<< HEAD
					MODEL_VP0,MODEL_RHO0,CXPB,MUTPB,TOURNSIZE)


def init_worker(base_seed):
    """
    Dá a cada processo do Pool um gerador de números aleatórios exclusivo.
    Isso é vital caso sua função `inversion_objective` use cálculos estocásticos.
    """
    global worker_rng
    
    pid = multiprocessing.current_process().pid
    sq = np.random.SeedSequence(entropy=base_seed, spawn_key=(pid,))
    worker_rng = np.random.default_rng(sq)

# ----------------------------------------------------------------------------------------------=

def configure_deap_ocean(vp_base, s_obs, seismic_config, dvp_min, dvp_max, num_layers, mutpb, tournsize, map_func,rng):
    """
    Configures the Genetic Algorithm for the Oceanographic Inversion of Delta Vp.
    Utilizes DEAP's native and bounded functions for maximum computational performance.

    Parameters
    ----------
    vp_base : array_like
        The initial background P-wave velocity model (e.g., from GLORYS climatology).
    s_obs : array_like
        The observed (or synthetic target) seismic trace data.
    seismic_config : dict
        Dictionary containing seismic forward modeling parameters (e.g., dt, z_max, dz, f_peak).
    dvp_min : float
        Lower physical bound for the velocity perturbation (Delta Vp) in m/s.
    dvp_max : float
        Upper physical bound for the velocity perturbation (Delta Vp) in m/s.
    num_layers : int
        Number of vertical layers in the 1D model (defines the number of genes).
    cxpb : float
        Probability of mating (crossover) between two individuals.
    mutpb : float
        Probability of mutation for each individual gene.
    tournsize : int
        Number of individuals participating in each tournament selection.
    map_func : callable, optional
        Mapping function for parallel execution (multiprocessing). Defaults to the built-in `map`.

    Returns
    -------
    toolbox : deap.base.Toolbox
        The fully configured DEAP toolbox ready for the evolutionary loop.
    """
    
    # 1. DNA Creation (Objective: Minimize Phase Misfit)
    if not hasattr(creator, "FitnessMin"):
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    if not hasattr(creator, "Individual"):
        creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()

    # 2. Parallelization (Multiprocessing mapping)
    toolbox.register("map", map_func)

    # 3. Initial Population Generation (Uniform Distribution)
    toolbox.register("attr_float", rng.uniform, dvp_min, dvp_max)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, n=num_layers)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # 4. Evaluation Function (Forward Modeling and NCC Calculation)
    # Note: Ensure the target function (e.g., 'phase_misfit_objective') is defined in English as well.
    toolbox.register("evaluate", phase_misfit_objective, 
                     base_model=vp_base, 
                     s_obs=s_obs, 
                     seismic_config=seismic_config)

    # --------------------
    # 5. Genetic Operators
    
    # Crossover: Simulated Binary Bounded (Weighted average respecting boundaries)
    toolbox.register("mate", tools.cxSimulatedBinaryBounded, 
                     low=dvp_min, up=dvp_max, eta=20.0)

    # Mutation: Polynomial Bounded (Local noise addition respecting boundaries)
    toolbox.register("mutate", tools.mutPolynomialBounded,low=dvp_min, up=dvp_max, eta=20.0, indpb=mutpb)

    # Selection: Tournament Strategy
    toolbox.register("select", tools.selTournament, tournsize=tournsize)
    
    return toolbox

# ----------------------------------------------------------------------------------------------
=======
					MAX_TOTAL,FREQ,MIN_FREQ,CXPB,MUTPB,TOURNSIZE)


def init_worker(base_seed):
    """
    Dá a cada processo do Pool um gerador de números aleatórios exclusivo.
    Isso é vital caso sua função `inversion_objective` use cálculos estocásticos.
    """
    global worker_rng
    
    pid = multiprocessing.current_process().pid
    sq = np.random.SeedSequence(entropy=base_seed, spawn_key=(pid,))
    worker_rng = np.random.default_rng(sq)
    
# ---------------------------------------------------------------------

def create_layers(min_thick_layer, max_thick_layer, max_total, max_layers, rng):
    '''
    This function generates a list of EXACTLY `max_layers` random layer thicknesses 
    that sum exactly to `max_total`. 
    
    Each layer thickness is dynamically bounded to ensure that all layers meet 
    both the `min_thick_layer` and `max_thick_layer` requirements.

    Parameters:
    -----------
    min_thick_layer : float
        Minimum thickness of an individual layer (m).
    max_thick_layer : float
        Maximum thickness of an individual layer (m).
    max_total : float
        Maximum total thickness of all layers (m).
    max_layers : int
        Exact number of layers of the model.
    rng : numpy.random.Generator
        A NumPy random number generator instance. 

    Returns:
    --------
    thick_lst : list of float
        A list of exactly `max_layers` layer thicknesses.
    '''

    values = []
    remaining = max_total
    
    for i in range(max_layers - 1):
        layers_left = max_layers - 1 - i
        
        lower_bound = max(min_thick_layer, remaining - (layers_left * max_thick_layer))
        
        upper_bound = min(max_thick_layer, remaining - (layers_left * min_thick_layer))
        
        choice = round(rng.uniform(lower_bound, upper_bound), 2)
        
        values.append(choice)
        remaining = round(remaining - choice, 2)
        
    values.append(round(remaining, 2))
    
    return values

# ------------------------------------------------------------------

def uniform(low_thick, up_thick, max_total, max_layers, low_vels, up_vels, rng):
    """
    Generates a synthetic 1D layered seismic velocity model for MASW inversion,
    where each layer has its own interpolated Vs bounds (top-to-bottom).

    Parameters
    ----------
    low_thick : float
        Minimum thickness allowed for a single layer.
    up_thick : float
        Maximum thickness allowed for a single layer.
    max_total : float
        Maximum total depth (cumulative thickness) of the profile.
    max_layers : int
        Maximum number of geological layers to generate.
    low_vels : list or tuple
        [min_surface_Vs, max_surface_Vs] — bounds for the shallowest layer.
    up_vels : list or tuple
        [min_deep_Vs, max_deep_Vs] — bounds for the deepest layer.
    rng : numpy.random.Generator
        Random number generator instance.

    Returns
    -------
    list of list
        [[thickness, Vs], ...] from top to bottom.
    """

    # -----------------
    # Layer thicknesses
    
    thickness_lst = create_layers(min_thick_layer=low_thick,max_thick_layer=up_thick,max_total=max_total,max_layers=max_layers,rng=rng)

    n = len(thickness_lst)

    # -------------------------
    # Per-layer velocity bounds
    # Linearly interpolated from surface → depth
    # t=0 → low_vels bounds
    # t=1 → up_vels bounds

    t = np.linspace(0, 1, n) # one scalar per layer

    low_bounds  = low_vels[0] + t * (up_vels[0] - low_vels[0])   # lower edge
    high_bounds = low_vels[1] + t * (up_vels[1] - low_vels[1])   # upper edge

    # -------------------------------
    # Sample each layer independently 
    # within its own [low, high] window
    vs_lst = [round(rng.uniform(lo, hi)) for lo, hi in zip(low_bounds, high_bounds)]

    # 4. Assemble model
    layer_model = [[t, v] for t, v in zip(thickness_lst, vs_lst)]

    return layer_model

# ----------------------------------------------------------------------------------------------

def inversion_objective(individual, true_disp, number_samples):
    """
    Objective function for inversion using DEAP with built-in physical constraints.

    This function evaluates the misfit between the experimental Rayleigh wave 
    dispersion data and the theoretical dispersion curve simulated from a given 
    earth model profile.

    Parameters
    ----------
    individual : list or array
        Estimated earth model profile used for optimization. Each layer is 
        expected to contain [thickness, Vs, Vp/Vs].
    true_disp : array
        Experimental Rayleigh wave phase velocity dispersion data.
    number_samples : int
        Number of frequency samples considered for misfit calculation.

    Returns
    -------
    tuple of float
        A single-element tuple containing the computed misfit value or a heavy penalty.
    
    Notes
    -----
    - Forward Modeling: The theoretical dispersion curve is obtained by creating a 
      velocity model from the Vs, then estimating Rayleigh phase velocities.
    - Misfit Formulation (RMSE):
      misfit = sqrt( sum((xdi - xci)^2) / nf )
      where:
        xdi: Experimental Rayleigh wave phase velocity.
        xci: Theoretical Rayleigh wave phase velocity.
        nf: Number of frequency samples.
    """

    PENALTY_OFFSET = 1e6

    # FORWARD MODELING & MISFIT CALCULATION
    try:
        simulated_velocity_model = create_velocity_model_from_profile_vs(individual)
            
        simulated_cpr = estimate_disp_from_velocity_model(
            vel_mol=simulated_velocity_model, 
            min_freq=MIN_FREQ, 
            max_freq=FREQ*3, 
            number_samples=number_samples)
            
        simulated_dispersion = simulated_cpr.velocity[::-1] * 1000

        misfit = np.sqrt(np.sum((true_disp - simulated_dispersion) ** 2) / number_samples)

        return misfit,

    except Exception as e:
        # Returns a high misfit if forward solver fails
        return PENALTY_OFFSET,

# ------------------------------------------------------------------
>>>>>>> 028cf37d1dc8a459a17f99af11a917a3fc5fb03d

def statistics_save(individual):
    """
    Retrieves the fitness value of an individual.

    This function returns the fitness value(s) of the given individual, 
    which is used for tracking statistics such as mean, standard deviation, 
    minimum, and maximum fitness during the optimization process.

    Parameters:
    -----------
    individual : object
        An individual solution with an assigned fitness value.

    Returns:
    --------
    fitness_value : tuple
        The fitness value(s) of the individual.
    """
    
    return individual.fitness.values

# ------------------------------------------------------------------

<<<<<<< HEAD
=======
def mutate_gaussian(ind,mutpb,rng):
    """
    Applies Gaussian mutation to the individual's layers.
    
    Parameters:
    -----------
    ind : list
        The individual consisting of layers, where each layer is a sublist 
        (e.g., [thickness, Vs, Vp/Vs]).
    mutpb : float
        Probability of mutating each value.
    rng : numpy.random.Generator
        A NumPy random number generator instance.
        
    Returns:
    --------
    tuple
        The mutated individual.
    """
    for i in range(len(ind)):  # Iterate over LAYERS
        for j in range(len(ind[i])):  # Iterate over PROPERTIES in the layer
            if rng.random() < mutpb:  # Mutation probability check
                value = ind[i][j]
                sigma = 0.1 * abs(value)  # Standard deviation as 10% of the current value
                ind[i][j] += round(rng.normal(0, sigma), 2)  # Apply Gaussian noise

    return ind,

# ------------------------------------------------------------------

def crossover_two_point(ind1, ind2, cxpb, rng):
    """
    Applies Two-Point Crossover by swapping whole layers as intact blocks.
    This preserves the physical relationship between a layer's thickness and velocities.
    """
    if rng.random() < cxpb:  # Check if crossover occurs for the whole individual
        num_layers = min(len(ind1), len(ind2))
        
        if num_layers > 1:
            point1, point2 = sorted(rng.choice(num_layers, size=2, replace=False))
            
            # Swap entire layers (thickness, vs) as a single unit
            for l in range(point1, point2 + 1):
                ind1[l], ind2[l] = ind2[l], ind1[l]
                
    return ind1, ind2

# ------------------------------------------------------------------

def configure_deap(map_func,lower_thick,upper_thick,max_total,max_layers,lower_vs,upper_vs,rng,estimated_disp,nb_samples):
    """
    Configures the DEAP toolbox for the evolutionary inversion of seismic velocity models.

    This function sets up the Genetic Algorithm (GA) environment. It defines the 
    minimization objective, initializes multiprocessing, maps the attribute generators 
    (for layer thicknesses, Vs, and Vp/Vs), and registers the evolutionary operators 
    (evaluation, crossover, mutation, and selection).

    Parameters
    ----------
    map_func : int
        Number of processes to use for parallel evaluation of individuals.
    lower_thick : float
        Lower bound for individual layer thickness (in meters).
    upper_thick : float
        Upper bound for individual layer thickness (in meters).
    max_total : float
        Maximum total thickness of all layers combined (in meters).
    max_layers : int
        Maximum number of layers allowed in the model. 
    lower_vs : float
        Lower bound for S-wave velocities (in m/s).
    upper_vs : float
        Upper bound for S-wave velocities (in m/s).
    rng : numpy.random.Generator
        A NumPy random number generator instance for reproducible stochastic operations.
    estimated_disp : array-like
        The true or target dispersion curve data used as the reference for the 
        inversion objective.
    nb_samples : int
        The number of samples (e.g., frequency or period points) in the dispersion curve.

    Returns
    -------
    toolbox : deap.base.Toolbox
        The fully configured DEAP toolbox, ready to run the evolutionary algorithm.

    Notes
    -----
    - This function relies on global variables `CXPB`, `MUTPB`, and `TOURNSIZE` 
      for crossover probability, mutation probability, and tournament size, respectively.
    - It assumes the existence of custom external functions: `uniform`, 
      `inversion_objective`, `crossover_two_point`, and `mutate_gaussian`.
    """

    # Fitness and Individual Creation:
    if not hasattr(creator, "FitnessMin"):
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    if not hasattr(creator, "Individual"):
        creator.create("Individual", list, fitness=creator.FitnessMin)

    # Toolbox Initialization:
    toolbox = base.Toolbox()

    # Using Multiple Processors 
    toolbox.register("map", map_func)
    
    # Attribute Generator:
    toolbox.register("model", uniform, low_thick=lower_thick, up_thick=upper_thick, max_total=max_total, max_layers=max_layers, low_vels=lower_vs, up_vels=upper_vs, rng=rng)

    # Individual and Population Initialization:
    toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.model)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # Evaluation Function:
    # L2 norm
    toolbox.register("evaluate", inversion_objective, true_disp=estimated_disp,number_samples=nb_samples)

    # Crossover Operation:
    toolbox.register("mate", crossover_two_point, cxpb=CXPB, rng=rng)

    # Mutation Operation:
    toolbox.register("mutate", mutate_gaussian, mutpb=MUTPB, rng=rng)
    
    # Selection Strategy:
    toolbox.register("select", tools.selTournament, tournsize=TOURNSIZE)
    
    return toolbox

# ------------------------------------------------------------------

def process_depths(thick_row):
    """
    Process depth list based on thickness values with special handling for MAX_TOTAL threshold.
    
    Parameters:
    -----------
    thickness_list : list or array-like
        List of layer thickness values from inversion results
    
    Returns:
    --------
    list
        Processed depth list with special handling for MAX_TOTAL threshold
    """
    # Calculate cumulative depths
    depths =  list(accumulate(thick_row[:-1])) + [2.0]

    return depths
# ------------------------------------------------------------------------------------

>>>>>>> 028cf37d1dc8a459a17f99af11a917a3fc5fb03d
def bootstrap_hof_uncertainty(df_input,rng, n_iterations=500, ci_percentiles=[2.5, 97.5]):
    """
    Perform bootstrap uncertainty analysis on Hall of Fame (HOF) solutions each 
    survey.
    
    Parameters:
    -----------
    df : dataframe
        Collection of best solutions from evolutionary algorithm (Hall of Fame)
    n_iterations : int, optional
        Number of bootstrap resamples (default: 1000)
    ci_percentiles : list, optional
        Percentiles for confidence intervals (default: [2.5, 97.5] for 95% CI)
    
    Returns:
    --------
    Dataframe
        Dataframe containing bootstrap statistics and confidence intervals
    """
  
    # Initialize storage for bootstrap statistics]
    survey_df = df_input['survey'].values[0]
    hof_vs = df_input['vs'].values
    hof_depth = df_input['depth'].values
    n_hof = len(hof_vs)

    bootstrap_vs_means = []
    bootstrap_depth_means = []

    # Bootstrap resampling process
    for i in range(n_iterations):
       
        # Resample with replacement
        sample = rng.choice(range(n_hof), size=n_hof, replace=True)

        bootstrap_vs = [hof_vs[sl] for sl in sample]
        bootstrap_depth = [hof_depth[sl] for sl in sample]
       
        # Calculate and store mean of resampled solutions
        bootstrap_vs_means.append(np.mean(bootstrap_vs, axis=0))
        bootstrap_depth_means.append(np.mean(bootstrap_depth, axis=0))
    
    # Calculate confidence intervals
    lower_vs, upper_vs = np.percentile(bootstrap_vs_means, ci_percentiles, axis=0)
    lower_depth, upper_depth = np.percentile(bootstrap_depth_means, ci_percentiles, axis=0)
    
    # Calculate statistics
    overall_mean_vs = np.mean(bootstrap_vs_means, axis=0)
    std_dev_vs = np.std(bootstrap_vs_means, axis=0)

    overall_mean_depth = np.mean(bootstrap_depth_means, axis=0)
    std_dev_depth = np.std(bootstrap_depth_means, axis=0)
    
    dic_bootstrap = {
        'survey': survey_df,
        'mean_vs': overall_mean_vs,
        'std_vs': std_dev_vs,
        'ci_lower_vs': lower_vs,
        'ci_upper_vs': upper_vs,
        'bootstrap_distribution_vs': bootstrap_vs_means,
        'mean_depth': overall_mean_depth,
        'std_depth': std_dev_depth,
        'ci_lower_depth': lower_depth,
        'ci_upper_depth': upper_depth,
        'bootstrap_distribution_depth': bootstrap_depth_means
    }

    return dic_bootstrap

# ------------------------------------------------------------------------------------
