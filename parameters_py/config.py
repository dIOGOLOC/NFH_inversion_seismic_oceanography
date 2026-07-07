"""
--------------------------------------------------------------------------------
         Module that parses global parameters from a configuration file
--------------------------------------------------------------------------------

Author: Diogo L.O.C. (locdiogo@gmail.com)


Last Date: 07/2026

Description:
Module that parses global parameters from a configuration file at first import,
to make them available to the other parts of the program.

More information in:
https://wiki.python.org/moin/ConfigParserExamples

Input:
Configuration file, wherein global paths and parameters are defined.

Outputs:
The module provides a parser for simple configuration files consisting of groups
of named values.

"""

import configparser
import os
import glob
import ast

def select_and_parse_config_file(basedir='.', ext='cnf', verbose=True):
    """
    Reads a configuration file and returns an instance of ConfigParser:
    First, looks for files in *basedir* with extension *ext*.
    Asks user to select a file if several files are found,
    and parses it using ConfigParser module.
    @rtype: L{ConfigParser.ConfigParser}
    """
    config_files = glob.glob(os.path.join(basedir, u'*.{}'.format(ext)))

    if not config_files:
        raise Exception("No configuration file found!")

    if len(config_files) == 1:
        # only one configuration file
        config_file = config_files[0]
    else:
        print("Select a configuration file:")
        for i, f in enumerate(config_files, start=1):
            print("{} - {}".format(i, f))
        res = int(input(''))
        config_file = config_files[res - 1]

    if verbose:
        print("Reading configuration file: {}".format(config_file))

    conf = configparser.ConfigParser(allow_no_value=True)
    conf.read(config_file)

    return conf

# ==========================
# parsing configuration file
# ==========================

config = select_and_parse_config_file(basedir='.', ext='cnf', verbose=True)

# -----
# paths
# -----

## ------------------
## Name of the folder

FOLDER_NAME = config.get('paths', 'FOLDER_NAME')

## -----------------------
## Directory of the output (Figures and Files)

FOLDER_OUTPUT = config.get('paths', 'FOLDER_OUTPUT')

# -----
# clima
# -----

## -------------------------------
## Global Ocean Physics Reanalysis - https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description
## This product includes daily and monthly mean files for temperature, salinity, currents, sea level, mixed layer depth 
## and ice parameters from the top to the bottom. 
## The global ocean output files are displayed on 50 standard levels. 

## ---------------------------------------------------
## File with the initial sound speed model from Glorys
## (Note: This should point to the saved numpy/netcdf file of the specific month, already interpolated)
MODEL_VP0 = config.get('clima', 'MODEL_VP0')

<<<<<<< HEAD
## -----------------------------------------------
## File with the initial density model from Glorys
## (Note: This should point to the saved numpy/netcdf file of the specific month, already interpolated)
MODEL_RHO0 = config.get('clima', 'MODEL_RHO0')


# -------
# sismica
# -------

## -----------------------------------
## Parameters for the Forward Modeling (1D Convolution)

## --------------------------------------------------
## Temporal sampling rate of the synthetic seismogram (e.g., 0.002 = 2 ms)
DT = config.getfloat('sismica', 'DT')

## ----------------------------------------------
## Maximum depth for the inversion grid in meters (defines the TWT window)
Z_MAX = config.getfloat('sismica', 'Z_MAX')

## ----------------------
## Grid spacing in meters (e.g., 20.0 m to match the 50 parameters over 1000 m)
DZ = config.getfloat('sismica', 'DZ')

## ---------------------------------------------------------------
## Central frequency of the source wavelet / low-pass filter in Hz
F_PICO = config.getfloat('sismica', 'F_PICO')
=======
## --------------
## Density ranges (g/cm³)
S_VELOCITY_RANGES = ast.literal_eval(config.get('model', 'S_VELOCITY_RANGES'))

## -----------
## VpVs ranges (dimensionless)
VPVS_RANGES = ast.literal_eval(config.get('model', 'VPVS_RANGES'))

## -------------------------
## S-wave velocity  basement (m/s²)
S_VELOCITY_BASEMENT = config.getfloat('model', 'S_VELOCITY_BASEMENT')

## -------------
## VpVs basement (dimensionless)
VPVS_BASEMENT = config.getfloat('model', 'VPVS_BASEMENT')

## -------------------------------------
## Maximum total thickness of all layers (m)
MAX_TOTAL = config.getfloat('model', 'MAX_TOTAL')

## -------------------------------
## Random seed for reproducibility
SEED = config.getint('model', 'SEED')

# ----------
# propagator
# ----------

## ---------------------------
## DEEPWAVE Elastic Propagator - https://ausargeo.com/deepwave/elastic
## The elastic propagator has three model parameters (the P and S wavespeeds and the Density)
## Deepwave internally uses lambda, mu ,and buoyancy (the two Lamé parameters and the reciprocal of density).

## -----
## MODEL

## Horizontal grid cell size in meters [default=2 cm]
DX = config.getfloat('propagator', 'DX')

## Vertical grid cell size in meters [default=2 cm]
DZ = config.getfloat('propagator', 'DZ')

## Total number of horizontal grid points [default=5000 points]
NX = config.getint('propagator', 'NX')

## Total number of vertical grid points [default=100 points]
NZ = config.getint('propagator', 'NZ')

## ---
## PML

## Original: top=1 (free surface), left/right/bottom = absorbing
## Deepwave pml_width = [top, bottom, left, right]
## Number of absorbing boundary points

PML_WIDTH = ast.literal_eval(config.get('propagator', 'PML_WIDTH'))

## -------
## SOURCES

## Total number of independent seismic shots to simulate
N_SHOTS = config.getint('propagator', 'N_SHOTS')

## Number of sources activated simultaneously per shot.
N_SHOTS_PER_SHOT = config.getint('propagator', 'N_SHOTS_PER_SHOT')

## Horizontal spacing between consecutive shot locations in grid points.
D_SOURCE = config.getint('propagator', 'D_SOURCE')

## Horizontal grid index for the first source location.
FIRST_SOURCE = config.getint('propagator', 'FIRST_SOURCE')

## Vertical grid index for the source (0 for surface).
SOURCE_DEPTH = config.getint('propagator', 'SOURCE_DEPTH')

## ---------
## RECEIVERS

## Total number of active receivers per shot.
N_RECEIVERS_PER_SHOT = config.getint('propagator', 'N_RECEIVERS_PER_SHOT')

## Horizontal spacing between adjacent receivers in grid points.
D_RECEIVERS = config.getint('propagator', 'D_RECEIVERS') 

## Horizontal grid index for the first receiver (sets the initial offset).
FIRST_RECEIVERS = config.getint('propagator', 'FIRST_RECEIVERS')

## Vertical grid index for the receiver array (0 for surface).
RECEIVERS_DEPTH = config.getint('propagator', 'RECEIVERS_DEPTH')

## -------
## WAVELET  

## Central frequency of the Ricker source wavelet in Hertz.
FREQ = config.getfloat('propagator', 'FREQ')

## Total number of discrete time steps to simulate and record.
NT = config.getint('propagator', 'NT')

## Time step duration in seconds, set to satisfy the CFL stability condition.
DT = config.getfloat('propagator', 'DT')

## The finite difference order of accuracy. Possible values are 2 and 8. [default=8].
ACCURACY = config.getint('propagator', 'ACCURACY')

## -------------------
## DISPERSION ANALYSIS  

# Minimum frequency threshold [Hz]. Frequencies below this are discarded to 
# prevent picking low-frequency noise artifacts. Default is 50.0.
MIN_FREQ = config.getfloat('propagator', 'MIN_FREQ')

# Standard deviation for the Gaussian filter used to smooth the image. 
# Default is 1.0.
SIGMA = config.getfloat('propagator', 'SIGMA')
>>>>>>> 028cf37d1dc8a459a17f99af11a917a3fc5fb03d

# ----
# gene
# ----

<<<<<<< HEAD
=======
## ----------------------------------------
## Minimum thickness of an individual layer (m).
MIN_THICK_LAYER = config.getfloat('gene', 'MIN_THICK_LAYER')

## ----------------------------------------
## Maximum thickness of an individual layer  (m)
MAX_THICK_LAYER = config.getfloat('gene', 'MAX_THICK_LAYER')

## -----------------------------
## Number of layers of the model (int)
MAX_LAYERS = config.getint('gene', 'MAX_LAYERS')

## ---------------------------------------
## First layer bounds for Vs values (m/s).
LOW_VELS = ast.literal_eval(config.get('gene', 'LOW_VELS'))

## --------------------------------------
## Last layer bounds for Vs values (m/s).
UP_VELS = ast.literal_eval(config.get('gene', 'UP_VELS'))

## ---------------
## First layer bound for Vp/Vs values.
FL_VPVS = config.getfloat('gene', 'FL_VPVS')

## ---------------
## Last layer bound for Vp/Vs values.
LL_VPVS = config.getfloat('gene', 'LL_VPVS')

>>>>>>> 028cf37d1dc8a459a17f99af11a917a3fc5fb03d
## ----------------------------------
## Probability of mutating each value (default=0.02)
MUTPB = config.getfloat('gene', 'MUTPB')

## ---------------------------------------
## The probability of performing crossover (default=0.7)
CXPB = config.getfloat('gene', 'CXPB')     

## ----------------------------------------------------------
## The number of individuals participating in each tournament
TOURNSIZE = config.getint('gene', 'TOURNSIZE')        

## ---------------------------------------
## Starting best solution estimation:
HOF_NUM = config.getint('gene', 'HOF_NUM')    

## ---------------------
## Number of generations:
NGEN = config.getint('gene', 'NGEN')

## ---------------
## Number of individuals in population:
POPULATION = config.getint('gene', 'POPULATION')

## --------------------
## Number of inversions:
N_INV = config.getint('gene', 'N_INV')

## ---------------
## MULTIPROCESSING
NUM_PROCESS = config.getint('gene', 'NUM_PROCESS')