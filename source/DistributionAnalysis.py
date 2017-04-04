"""
@author Tomas Lazauskas, 2017

A module to analyse distribution functions of different stoichiometries

"""

import numpy as np

import Constants
import Utilities
from Utilities import log

def calc_permutation(m, mm, _accuracy):
  """
  Evaluates the permutation expression
  
  """
  
  cat_sites = 108
  an_sites = 216
  
  up_fac = cat_sites - 2*m
  down_fac = cat_sites - 2*mm
  
  value = Utilities.factorial_devision(up_fac, down_fac, _accuracy)
  
  value *= Utilities.factorial_devision(2*m, 2*mm, _accuracy)
  
  up_fac = an_sites - m
  down_fac = an_sites - mm
  
  value *= Utilities.factorial_devision(up_fac, down_fac, _accuracy)
  
  value *= Utilities.factorial_devision(m, mm, _accuracy)
  
  return value

def distribution_analysis(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, 
                        experiment_cnts, _accuracy, options):
  """
  """
  
  success = True
  error = ""
  
  # Preparing Wm probabilities    
  Wm_array = prepare_Wm(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, 
                        experiment_cnts, _accuracy, options)
  
  # Plotting the Wm probabilities
  
  return success, error

def prepare_Wm(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, experiment_cnts, 
               _accuracy, options):
  """
  Evaluates Wm with respect to temperature and chemical potential
  
  """
  
  log(__name__, "Preparing Wm", options.verbose, indent=3)
  
  temp_len = len(temperatures)
  chem_pot_len = len(chem_pot_range)
  chem_pot_multi_len = len(chem_pot_multi)
  
  Wm_array = np.empty([chem_pot_multi_len, temp_len, chem_pot_len], _accuracy)
  
  for m_idx in range(chem_pot_multi_len):
    m_value = chem_pot_multi[m_idx]
    Emin_m = min_energies[m_idx]
    
    # Estimating the top part of Wm
    sum_top = delta_E_sums[m_idx]
    
    for t_i in range(temp_len):
      kT = np.longdouble(Constants.kB * temperatures[t_i])
     
      mu_i = 0
      for mu in chem_pot_range:    
        
        Wm_value = _accuracy(0.0)
        
        # Estimating the bottom part of Wm
        sum_bottom = _accuracy(0.0)
        
        for mm_index in range(chem_pot_multi_len):
          mm_value = chem_pot_multi[mm_index]
          Emin_mm = min_energies[mm_index]
          
          # exponential term
          exp_expr_pow = _accuracy(((Emin_m - Emin_mm) + mu*(m_value - mm_value))/kT)
          exp_expr = _accuracy(np.exp(exp_expr_pow))
          
          # Nm/Nmm
          attempts_cnt = experiment_cnts[m_idx] / experiment_cnts[mm_index] 
          
          # Pmm/Pm
          permutation = _accuracy(calc_permutation(m_value, mm_value, _accuracy))

          sum_bottom += exp_expr * attempts_cnt * permutation * delta_E_sums[mm_index]
        
        Wm_value = sum_top / sum_bottom
        
        Wm_array[m_idx, t_i, mu_i] = Wm_value
        
        mu_i += 1
      
  return Wm_array