"""
@author Tomas Lazauskas, 2017-2018

A module to analyse distribution functions of different stoichiometries

"""

import numpy as np
import sys

import Constants
import IO
import Graphs
import Utilities
from Utilities import log

# TODO: needs to made universal
def calc_permutation(m, mm, _accuracy):
  """
  Evaluates the permutation expression
  
  """
    
  value = 1
  
  return value

def average_analysis(temperatures, chem_pot_range, prop_array, prop_name, Wm_array, _accuracy, options, temp_depend=False):
  """
  Performs average analysis on a selected property
  
  """
  
  # Calculates the average values
  success, error, avg_array = calc_average_value(temperatures, chem_pot_range, prop_array, prop_name, 
                                                 Wm_array, _accuracy, options, temp_depend=temp_depend)
  
  if not success:
    return success, error
  
  # Plots the average values
  Graphs.avg_values(temperatures, chem_pot_range, avg_array, prop_name, _accuracy, options)
  
  return success, error

def calc_average_value(temperatures, chem_pot_range, prop_array, prop_name, Wm_array, _accuracy, options, temp_depend=False):
  """
  Calculates average value of a system's property
  
  """
  
  success = True
  error = ""
  
  log(__name__, "Calculating an average value of: %s" % (prop_name), options.verbose, indent=3)
    
  temp_len = len(temperatures)
  chem_pot_len = len(chem_pot_range)
  
  # is the value temperature dependent
  if not temp_depend:
    prop_len = len(prop_array)
  else:
    prop_len = len(prop_array[0])
  
  avg_array = np.zeros([temp_len, chem_pot_len], _accuracy)
    
  # for each temperature:
  for t_index in range(temp_len):
    temperature = temperatures[t_index]
    
    for mu_index in range(chem_pot_len):   
      
      prop_avg = _accuracy(0.0)
      
      for prop_index in range(prop_len):
        wm_value = Wm_array[t_index, mu_index, prop_index]
        
        # is the value temperature dependent
        if not temp_depend:
          prop_value = prop_array[prop_index]
        else:
          prop_value = prop_array[t_index][prop_index]
        
        prop_avg += wm_value * prop_value
      
      avg_array[t_index, mu_index] = prop_avg

  return success, error, avg_array

def distribution_analysis(chem_pot_multi, names, temperatures, chem_pot_range, min_energies, delta_E_sums, 
                        experiment_cnts, permutations, omega_c_arr, _accuracy, options):
  """
  Performs the distribution analysis: evaluates Wm and plots it against m and mu.
  
  """
  
  log(__name__, "Distribution analysis", options.verbose, indent=2)
  
  success = True
  error = ""
  
  # Preparing Wm probabilities    
  Wm_array = prepare_Wm(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, 
                        experiment_cnts, permutations, _accuracy, options)
  
  # Writing Wm into a file
  success, error = IO.write_Wm(temperatures, chem_pot_range, chem_pot_multi, Wm_array)
  
  if not success:
    return success, error, Wm_array
  
  # Plotting the Wm probabilities 3D plots
  Graphs.wm_contour(temperatures, names, chem_pot_range, chem_pot_multi, Wm_array, _accuracy, options)
  
  # Performing analysis with respect to the distribution function (average m is the standard analysis)
  # average m
  average_analysis(temperatures, chem_pot_range, chem_pot_multi, "m", Wm_array, _accuracy, options, temp_depend=False)
  
  # average gamma
  average_analysis(temperatures, chem_pot_range, omega_c_arr, "\gamma^{c}", Wm_array, _accuracy, options, temp_depend=True)
  
  return success, error, Wm_array

def distribution_analysis_canonical(temperatures, experiment_cnts, shifted_energies, delta_E_sums, _accuracy, options):
  """
  Evaluates Wm with respect to temperature
  
  """
  
  log(__name__, "Distribution analysis (canonical)", options.verbose, indent=2)
  
  success = True
  error = ""
  
  temp_len = len(temperatures)
  comp_len = len(experiment_cnts)
  max_experiment_cnt = np.int32(np.max(experiment_cnts))
  
  Wm_array = np.zeros([temp_len, comp_len, max_experiment_cnt], _accuracy)
  
  # for each temperature
  for t_i in range(temp_len):
    kT = np.longdouble(Constants.kB * temperatures[t_i])
    
    # for each discrete composition
    for comp_i in range(comp_len):
      
      # for each experiment
      for exp_i in range(experiment_cnts[comp_i]):
                
        exp_expr_pow = -1.0 * shifted_energies[comp_i][exp_i] / kT
        
        Wm_array[t_i, comp_i, exp_i] = np.exp(exp_expr_pow) / delta_E_sums[comp_i, t_i]
  
  return success, error, Wm_array

def prepare_Wm(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, experiment_cnts, 
               permutations, _accuracy, options):
  """
  Evaluates Wm with respect to temperature and chemical potential
  
  """
  
  log(__name__, "Preparing Wm", options.verbose, indent=3)
  
  temp_len = len(temperatures)
  chem_pot_len = len(chem_pot_range)
  chem_pot_multi_len = len(chem_pot_multi)
  
  Wm_array = np.zeros([temp_len, chem_pot_len, chem_pot_multi_len], _accuracy)
  
  for t_i in range(temp_len):
    kT = np.longdouble(Constants.kB * temperatures[t_i])
  
    for mu_i in range(chem_pot_len):    
      mu = chem_pot_range[mu_i]
      
      Wm_sum = _accuracy(0.0)
      for m_index in range(chem_pot_multi_len):
        
        m_value = chem_pot_multi[m_index]
        Emin_m = min_energies[m_index]
        
        # Estimating the top part of Wm
        sum_top = delta_E_sums[m_index][t_i]
        
        # Estimating the bottom part of Wm
        sum_bottom = _accuracy(0.0)
                
        breakLoops = False
  
        for mm_index in range(chem_pot_multi_len):

          mm_value = chem_pot_multi[mm_index]
          Emin_mm = min_energies[mm_index]
          
          # exponential terms
          
          temp_expr1 = (Emin_m - Emin_mm)
          #temp_expr1 = np.abs((Emin_m - Emin_mm))
          
          temp_expr2 = mu*(m_value - mm_value)
          #temp_expr2 = np.abs(mu*(m_value - mm_value))
          
          exp_expr_pow = _accuracy((temp_expr1 - temp_expr2)/kT)
                      
          if (exp_expr_pow > Constants.BIGEXPO):            
            exp_expr = Utilities.approximate_exp(exp_expr_pow, _accuracy)
            
          else:
            exp_expr = _accuracy(np.exp(exp_expr_pow))
          
          # Nm/Nmm
          attempts_cnt = experiment_cnts[m_index] / experiment_cnts[mm_index] 
          
          # Pmm/Pm
          if options.permCalc:
            permutation = _accuracy(calc_permutation(m_value, mm_value, _accuracy))
          else:
            permutation = _accuracy(_accuracy(permutations[mm_index])/_accuracy(permutations[m_index]))
            
          sum_bottom += (((attempts_cnt * delta_E_sums[mm_index][t_i]) * permutation) * exp_expr)
          
        if not breakLoops:
           
          if sum_bottom == 0.0:
            Wm_value = 1.1
          else:
            Wm_value = sum_top / sum_bottom
            
          Wm_sum += Wm_value
        else:
          Wm_value = 0
        
        Wm_array[t_i, mu_i, m_index] = Wm_value
       
      # Normalising
      if Wm_sum != 0.0:
        Wm_array[t_i, mu_i, :] /= Wm_sum
  
  return Wm_array
