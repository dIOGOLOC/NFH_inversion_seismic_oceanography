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
MES_ALVO = config.get("clima", "MES_ALVO")
N_EOFS = config.getint("clima", "N_EOFS")
 
# ---------------------------------------------------------------------
# [sismica]
# ---------------------------------------------------------------------

DT = config.getfloat("sismica", "DT")
Z_MAX = config.getfloat("sismica", "Z_MAX")
DZ = config.getfloat("sismica", "DZ")
F_PICO = config.getfloat("sismica", "F_PICO")
WAVELET_TYPE = config.get("sismica", "WAVELET_TYPE")
WAVELET_LENGTH = config.getfloat("sismica", "WAVELET_LENGTH")
 
# ---------------------------------------------------------------------
# [gene]
# ---------------------------------------------------------------------

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