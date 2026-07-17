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

def select_and_parse_config_file(basedir=None, ext="cnf", verbose=True):
    """
    Finds a single configuration file with extension `ext` and parses it.
 
    If `basedir` is None (default), searches the repository root -- i.e.
    the parent directory of this package (ocean_inversion/../) -- so the
    search does NOT depend on the current working directory.
    """
    if basedir is None:
        # Repo root = parent directory of this package (ocean_inversion/..)
        basedir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
 
    config_files = sorted(glob.glob(os.path.join(basedir, f"*.{ext}")))
 
    if not config_files:
        raise FileNotFoundError(
            f"No '*.{ext}' configuration file found in '{os.path.abspath(basedir)}'."
        )
 
    if len(config_files) == 1:
        config_file = config_files[0]
    else:
        print("Several configuration files were found:")
        for i, f in enumerate(config_files):
            print(f"  [{i}] {f}")
        idx = int(input("Select the config file number to use: "))
        config_file = config_files[idx]
 
    if verbose:
        print(f"[ocean_inversion.config] Reading configuration file: {config_file}")
 
    conf = configparser.ConfigParser(inline_comment_prefixes=("#",))
    conf.read(config_file)
    return conf
 
 
# basedir=None -> resolved automatically relative to this file's location,
# regardless of where Python/Jupyter was launched from.
config = select_and_parse_config_file(ext="cnf", verbose=True)
 

# ---------------------------------------------------------------------
# [paths]
# ---------------------------------------------------------------------

FOLDER_NAME = config.get("paths", "FOLDER_NAME")
FOLDER_OUTPUT = config.get("paths", "FOLDER_OUTPUT")
OBSERVED_DATA = config.get("paths", "OBSERVED_DATA")
 
# ---------------------------------------------------------------------
# [clima]
# ---------------------------------------------------------------------

MODEL_VP0 = config.get("clima", "MODEL_VP0")
MODEL_RHO0 = config.get("clima", "MODEL_RHO0")
MODEL_EOFS = config.get("clima", "MODEL_EOFS")
MODEL_COEF_HIST = config.get("clima", "MODEL_COEF_HIST")
MONTH_TARGET = config.get("clima", "MONTH_TARGET")
N_EOFS = config.getint("clima", "N_EOFS")
 
# ---------------------------------------------------------------------
# [sismica]
# ---------------------------------------------------------------------

DT = config.getfloat("sismica", "DT")
NT = config.getint("sismica", "NT")
Z_MAX = config.getfloat("sismica", "Z_MAX")
DZ = config.getfloat("sismica", "DZ")
F_PEAK = config.getfloat("sismica", "F_PEAK")
ADD_NOISE = config.getboolean("sismica", "ADD_NOISE")
PERCENTAGE_NOISE = config.getfloat("sismica", "PERCENTAGE_NOISE")

# ---------------------------------------------------------------------
# [gene]
# ---------------------------------------------------------------------

SEED = config.getint("gene", "SEED")
MUTPB = config.getfloat("gene", "MUTPB")
CXPB = config.getfloat("gene", "CXPB")
TOURNSIZE = config.getint("gene", "TOURNSIZE")
HOF_NUM = config.getint("gene", "HOF_NUM")
NGEN = config.getint("gene", "NGEN")
POPULATION = config.getint("gene", "POPULATION")
N_INV = config.getint("gene", "N_INV")
NUM_PROCESS = config.getint("gene", "NUM_PROCESS")
MARGEM_COEF = config.getfloat("gene", "MARGEM_COEF")
ETA_GENE = config.getfloat("gene", "ETA_GENE")