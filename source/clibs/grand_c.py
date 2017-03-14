# -*- coding: utf-8 -*-

"""
Wrapper to calculations.c

"""
import os

from ctypes import CDLL, c_double, c_int, POINTER

from .numpy_utils import CPtrToDouble, CPtrToInt

################################################################################

# load lib
_lib = CDLL(os.path.join(os.path.dirname(__file__), "_grand_c.so"))

################################################################################

# prototypes
_lib.calculate_distribution_functions.restype = c_int
_lib.calculate_distribution_functions.argtypes = [c_double, c_double, 
                                                  POINTER(c_double), 
                                                  POINTER(c_double), 
                                                  c_int, POINTER(c_double), 
                                                  POINTER(c_double)]

_lib.calculate_omega.restype = c_double
_lib.calculate_omega.argtypes = [c_double, c_double, c_double, c_double, 
                                 POINTER(c_double), 
                                 POINTER(c_double), 
                                 c_int, POINTER(c_double)]

# interfaces
def calculate_distribution_functions(temperature, mu, energies_mins, energies_diffs, 
                                     stoich_cnt, stoichiometries, wms):
  """
  Calculates distribution functions for all the stoichiometries.
  
  """
    
  return _lib.calculate_distribution_functions(temperature, mu, 
                                               CPtrToDouble(energies_mins),
                                               CPtrToDouble(energies_diffs), 
                                               stoich_cnt, CPtrToDouble(stoichiometries), 
                                               CPtrToDouble(wms))

def calculate_gamma(temperature, mu, energies_GM, area, energies_mins, energies_diffs, 
                                     stoich_cnt, stoichiometries):
  """
  Calculates omega with respect to temperature and mu.
  
  """
    
  return _lib.calculate_gamma(temperature, mu, energies_GM, area, 
                              CPtrToDouble(energies_mins),
                              CPtrToDouble(energies_diffs), 
                              stoich_cnt, CPtrToDouble(stoichiometries))