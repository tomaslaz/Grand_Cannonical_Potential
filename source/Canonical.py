#!/usr/bin/python
"""
@author Tomas Lazauskas, 2017

Canonical ensemble module.

"""

import DistributionAnalysis
import OmegaAnalysis
from Utilities import log
  
def perform_canonical_analysis(chem_pot_multi, names, permutations, temperatures, min_energies, 
                               delta_E_sums, experiment_cnts, shifted_energies, _accuracy, options):
  """
  The main method of the grand canonical analysis
  
  """
  
  success, error, Wm_array = DistributionAnalysis.distribution_analysis_canonical(temperatures, experiment_cnts, shifted_energies, 
                                                                                  delta_E_sums, _accuracy, options)
   
  if not success:
    return success, error
  
  # perform omega analysis
  success, error, = OmegaAnalysis.c_gamma_analysis(chem_pot_multi, names, temperatures, min_energies, 
                                                   delta_E_sums, experiment_cnts, permutations, _accuracy, options)
  
  return success, error