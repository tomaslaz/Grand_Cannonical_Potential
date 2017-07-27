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

def c_calc_gamma(temperatures, min_energies, delta_E_sums, experiment_cnts, permutations, _accuracy, options):
  """
  Calculates omega values with respect to temperature (canonical analysis)
  
  """
  
  success = True
  error = ""
  
  log(__name__, "Calculating Omega (canonical)", options.verbose, indent=3)
  
  temp_len = len(temperatures)
  comp_len = len(min_energies)
  
  omega_arr = np.zeros([temp_len, comp_len], _accuracy)
  
  # for each temperature:
  for t_index in range(temp_len):
    temperature = temperatures[t_index]
    
    kT = np.longdouble(Constants.kB * temperature)
    # for each composition
    
    for c_index in range(comp_len):
      # calculating Z^c_m
      
      # Pm
      if options.permCalc:
        Pm = _accuracy(calc_permutation(m_value, mm_value, _accuracy))
      else:
        Pm = _accuracy(permutations[c_index])         
      
      # Nm
      Nm = experiment_cnts[c_index]
      
      Z_cm = np.exp(-1.0*(min_energies[c_index]) / kT) * (Pm/Nm) * delta_E_sums[c_index][t_index]
      
      omega_value = - kT * np.log(Z_cm)
      
      omega_arr[t_index, c_index] = omega_value / Constants.gamma_c_m_coef
      
  return success, error, omega_arr

def c_gamma_analysis(chem_pot_multi, names, temperatures, min_energies, delta_E_sums, experiment_cnts, permutations, 
                     _accuracy, options):
  """
  Performs the canonical analysis.
  
  """
  
  success = True
  error = ""
  
  log(__name__, "Omega analysis (canonical)", options.verbose, indent=2)
  
  # calculates omega values
  success, error, omega_arr = c_calc_gamma(temperatures, min_energies, delta_E_sums, experiment_cnts, permutations, 
                                           _accuracy, options)
  
  if not success:
    return success, error
  
  # plot omega with respect to temperature
  Graphs.c_omega(chem_pot_multi, names, temperatures, omega_arr, _accuracy, options)

  return success, error

def g_c_calc_omega(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, experiment_cnts, 
                   permutations, _accuracy, options):
  """
  Calculates omega values with respect to temperature and chemical potential (grand canonical analysis)
  
  """
    
  success = True
  error = ""
  
  log(__name__, "Calculating Omega (grand canonical)", options.verbose, indent=3)
  
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
      
      # for each composition
      for m_index in range(chem_pot_multi_len):
        m_value = chem_pot_multi[m_index]
        
        min_energy = min_energies[m_index]
        
        exp_expr = _accuracy(-1.0*(min_energy - global_min_energy + m_value*mu_value) / (kT))
        
        sum2 += np.exp(exp_expr) * (permutations[m_index]/experiment_cnts[m_index]) * delta_E_sums[m_index][t_index]
            
      omega_value = -kT * (-global_min_energy + np.log(sum2))
              
      omega_arr[t_index, mu_index] = omega_value
  
  return success, error, omega_arr

def g_c_omega_analysis(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, experiment_cnts, 
                   permutations, _accuracy, options):
  """
  Performs the grand potential analysis.
  
  """
  
  success = True
  error = ""
  
  log(__name__, "Omega analysis (grand canonical)", options.verbose, indent=2)
  
  # calculates omega values
  success, error, omega_arr = g_c_calc_omega(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, 
                                         experiment_cnts, permutations, _accuracy, options)
  
  if not success:
    return success, error
  
  # plot omega over mu
  Graphs.c_g_omega_over_mu(temperatures, chem_pot_range, omega_arr, _accuracy, options)

  return success, error