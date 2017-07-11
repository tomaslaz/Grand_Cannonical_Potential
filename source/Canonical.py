#!/usr/bin/python
"""
@author Tomas Lazauskas, 2017

Canonical ensemble module.

"""

import OmegaAnalysis
from Utilities import log
  
def perform_canonical_analysis(chem_pot_multi, names, permutations, temperatures, min_energies, 
                               delta_E_sums, experiment_cnts, _accuracy, options):
  """
  The main method of the grand canonical analysis
  
  """
  
  success = True
  error = ""
  
  # perform omega analysis
  success, error, = OmegaAnalysis.c_omega_analysis(chem_pot_multi, names, temperatures, min_energies, 
                                                   delta_E_sums, experiment_cnts, permutations, _accuracy, options)
  
  return success, error