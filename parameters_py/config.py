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

# ----
# gene
# ----

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