"""
@author Tomas Lazauskas, 2017

A module to calculate and analyse Omega (the grand potential)

"""

import numpy as np
import sys

import Constants
import IO
import Graphs
import Utilities
from Utilities import log

def calc_omega(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, experiment_cnts, 
                   permutations, _accuracy, options):
  """
  Calculates omega values with respect to temperature and chemical potential
  
  """
    
  success = True
  error = ""
  
  log(__name__, "Calculating Omega", options.verbose, indent=3)
  
  temp_len = len(temperatures)
  chem_pot_len = len(chem_pot_range)
  chem_pot_multi_len = len(chem_pot_multi)
  
  global_min_energy = np.min(min_energies)
  
  omega_arr = np.zeros([temp_len, chem_pot_len], _accuracy)
  
  # for each temperature:
  for t_index in range(temp_len):
    temperature = temperatures[t_index]
    
    kT = np.longdouble(Constants.kB * temperature)
    
    # for each chemical potential value
    for mu_index in range(chem_pot_len):
      mu_value = chem_pot_range[mu_index]
      
      sum2 = _accuracy(0.0)
      
      # for each stoichiometry
      for m_index in range(chem_pot_multi_len):
        m_value = chem_pot_multi[m_index]
        
        min_energy = min_energies[m_index]
        
        sum2 += np.exp(-1.0*(min_energy - global_min_energy + m_value*mu_value) / (kT)) * delta_E_sums[m_index][t_index]
            
      omega_value = global_min_energy - kT * np.log(sum2)
              
      omega_arr[t_index, mu_index] = omega_value
  
  return success, error, omega_arr

def omega_analysis(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, experiment_cnts, 
                   permutations, _accuracy, options):
  """
  Performs the grand potential analysis.
  
  """
  
  success = True
  error = ""
  
  log(__name__, "Omega analysis", options.verbose, indent=2)
  
  # calculates omega values
  success, error, omega_arr = calc_omega(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, 
                                         experiment_cnts, permutations, _accuracy, options)
  
  if not success:
    return success, error
  
  Graphs.omega_over_mu(temperatures, chem_pot_range, omega_arr, _accuracy, options)

  return success, error