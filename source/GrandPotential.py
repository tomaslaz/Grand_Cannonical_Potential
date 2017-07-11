"""
@author Tomas Lazauskas, 2016-2017

Grand canonical ensemble module.
"""

import copy
import math
import numpy as np
import sys

import Constants
import Data
import DistributionAnalysis
import OmegaAnalysis

from Utilities import log
  
def perform_grand_canonical_analysis(permutations, chem_pot_multi, temperatures, chem_pot_range, 
                                     min_energies, delta_E_sums, experiment_cnts, _accuracy, options):
  """
  The main method of the grand canonical analysis
  
  """
    
  success, error, Wm_array = DistributionAnalysis.distribution_analysis(chem_pot_multi, temperatures, chem_pot_range, 
                                         min_energies, delta_E_sums, experiment_cnts, permutations, _accuracy, options)
  
  if not success:
    return success, error
  
  # perform omega analysis
  success, error = OmegaAnalysis.g_c_omega_analysis(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, 
                               experiment_cnts, permutations, _accuracy, options)
  
  return success, error
