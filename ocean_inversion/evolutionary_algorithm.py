# Functions 

## Dispersion curve estimative 

import numpy as np
import pandas as pd
from itertools import accumulate
import multiprocessing


from deap import base, creator, tools, algorithms

from ocean_inversion.modeling import phase_misfit_objective

from ocean_inversion.config import (
					MODEL_VP0,MODEL_RHO0,CXPB,MUTPB,TOURNSIZE)


def init_worker(base_seed):
    """
    Initialize a worker process with its own independent random number generator.

    Intended as the ``initializer`` for a multiprocessing.Pool: each worker
    process gets a dedicated NumPy Generator, seeded deterministically from
    a shared base seed combined with the worker's PID. This avoids
    correlated or duplicate random streams across parallel processes while
    keeping results reproducible for a given base_seed.

    Parameters
    ----------
    base_seed : int
        Base seed used to derive independent, reproducible seed sequences
        for each worker process.

    Notes
    -----
    Sets the module-level global variable `worker_rng`, which each worker
    process can then use for all subsequent random sampling.
    """
    global worker_rng
    
    pid = multiprocessing.current_process().pid
    sq = np.random.SeedSequence(entropy=base_seed, spawn_key=(pid,))
    worker_rng = np.random.default_rng(sq)


def create_individual_eof(mins, maxs, rng):
    """
    Generate an individual whose genes (EOF coefficients) respect their
    own bounds.

    Each gene is independently sampled from a uniform distribution within
    its corresponding [min, max] range.

    Parameters
    ----------
    mins : sequence of float
        Lower bound for each gene (EOF coefficient).
    maxs : sequence of float
        Upper bound for each gene (EOF coefficient), same length as mins.
    rng : numpy.random.Generator
        Random number generator used to draw the uniform samples.

    Returns
    -------
    list of float
        A new individual, with one value per gene, each drawn uniformly
        from its respective [mins[i], maxs[i]] range.
    """
    return [rng.uniform(mn, mx) for mn, mx in zip(mins, maxs)]

# ----------------------------------------------------------------------------------------------

def configure_deap_ocean(vp_base, eof_basis, coef_mins, coef_maxs, s_obs, f_peak, dt, dz, nt,z_max, rho_model, eta, mutpb, tournsize, map_func, rng):
    """
    Configures the Genetic Algorithm for the Oceanographic Inversion using EOFs.
    Utilizes DEAP's native and bounded functions for maximum computational performance.

    Parameters
    ----------
    map_func : int
        Number of processes to use for parallel evaluation of individuals.
    vp_base : array_like
        The initial background P-wave velocity model (MODEL_VP0).
    eof_basis : array_like
        The retained EOF basis matrix (MODEL_EOFS) of shape (K_reg, N_EOFS).
    coef_mins : list or array_like
        Lower physical bounds for each EOF coefficient.
    coef_maxs : list or array_like
        Upper physical bounds for each EOF coefficient.
    s_obs : array_like
        The observed (or synthetic target) seismic trace data.
    f_peak : float
        Central frequency of the source wavelet / low-pass filter in Hz
    dt : float
        Temporal sampling rate of the synthetic seismogram
    dz : float
        Temporal sampling rate of the synthetic seismogram
    nt : float
        Temporal sampling rate of the synthetic seismogram
    z_max : int
        Maximum depth for the inversion grid in meters (defines the TWT window).
    mutpb : float
        Probability of mutation for each individual gene.
    tournsize : int
        Number of individuals participating in each tournament selection.
    map_func : callable
        Mapping function for parallel execution (multiprocessing).
    rng : numpy.random.Generator
        Random number generator.

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

    # 3. Initial Population Generation (Specific limits for each EOF)
    toolbox.register("attr_genes", create_individual_eof, coef_mins, coef_maxs, rng)
    toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.attr_genes)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # 4. Evaluation Function
    toolbox.register("evaluate", phase_misfit_objective, 
                     base_model=vp_base,
                     eof_basis=eof_basis,
                     s_obs=s_obs, 
                     f_peak=f_peak, 
                     dt=dt, 
                     dz=dz, 
                     nt=nt,
                     z_max=z_max,  
                     rho_model=rho_model)

    # --------------------
    # 5. Genetic Operators
    
    # Crossover: Simulated Binary Bounded aceita listas para low e up
    toolbox.register("mate", tools.cxSimulatedBinaryBounded, 
                     low=coef_mins.tolist(), up=coef_maxs.tolist(), eta=eta)

    # Mutation: Polynomial Bounded aceita listas para low e up
    toolbox.register("mutate", tools.mutPolynomialBounded, 
                     low=coef_mins.tolist(), up=coef_maxs.tolist(), eta=eta, indpb=0.1)

    # Selection: Tournament Strategy
    toolbox.register("select", tools.selTournament, tournsize=tournsize)
    
    return toolbox

# ----------------------------------------------------------------------------------------------

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