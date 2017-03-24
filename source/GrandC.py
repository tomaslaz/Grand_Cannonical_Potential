"""
@author Tomas Lazauskas, 2016-2017

Grand canonical ensemble module.
"""

import copy
import math
import numpy as np
import sys

import Constants
import Utilities
from Utilities import log

def perform_grand_canonical_analysis(input_data_array, permutations, _accuracy, options):
  """
  The main method of the grand canonical analysis
  
  """
  
  success = True
  error = ""
  
  energy_list = []
  name_list = []
  total_combi_list = []
  
  # reading in the temperatures
  success, error, temperatures = Utilities.get_list_of_temps(options.temps)
  
  if not success:
    return success, error
  
  # reading in the chemical potential
  if options.urange is not None:
    success, error, chem_pot_range = Utilities.get_chem_pot_range(options.urange)
    
    print "chem_pot_range: ", chem_pot_range
    
    if not success:
      return success, error
    
  sys.exit()
  
  log(__name__, "Analysis will be performed at (K): %s" % (options.temps), options.verbose, indent=2)
  
  # first of all lets prepare the data for the calculations.
  log(__name__, "Preparing the data", options.verbose, indent=2)
  
  energies, min_energies, shifted_energies = prepare_energies(input_data_array, _accuracy, options)
  
  Z_c_m_sums = prepare_Z_c_m(shifted_energies, temperatures, permutations, _accuracy, options)
  
  print Z_c_m_sums
  
  return success, error

def prepare_energies(input_data_array, _accuracy, options):
  """
  Prepares energy arrays for the calculations
  
  """
  
  log(__name__, "Preparing the energies", options.verbose, indent=3)
  
  cases_cnt = input_data_array.shape[0]
  
  min_energies = np.empty(cases_cnt, dtype=_accuracy)
  energies = []
  shifted_energies = []
  
  for i in range(cases_cnt):
    # removes all zero values from the beginning and end of the array
    case_energies = np.trim_zeros(copy.deepcopy(input_data_array[i]))
    energies.append(case_energies)
    
    # it might be the first value but it might not be, thus we getting the minimum
    min_energies[i] = _accuracy(np.min(case_energies))
    
    diff_energy = case_energies - min_energies[i]
    shifted_energies.append(diff_energy)
      
  return energies, min_energies, shifted_energies

def prepare_Z_c_m(shifted_energies, temperatures, permutations, _accuracy, options):
  """
  Prepares the sum of the second part of the Z_c_m equation
  
  """
    
  energies_cnt = len(shifted_energies)
  temp_cnt = len(temperatures)
  
  Z_c_m_sums = np.empty([energies_cnt, temp_cnt], dtype=_accuracy)
    
  log(__name__, "Preparing Z_c_m_sums to speed up the calculations", options.verbose, indent=3)
  
  for t_i in range(temp_cnt):
    temperature = temperatures[t_i]
    
    kT = _accuracy(Constants.kB * temperature)
    
    for e_i in range(energies_cnt):
      
      case_energies = shifted_energies[e_i]
      cnt_case_energies = len(case_energies)
      
      sum_Ediff = _accuracy(0.0)
      
      # summing in the reverse order in order to account for small values
      for i in reversed(range(cnt_case_energies)): 
        sum_Ediff += _accuracy(math.exp(-1.0*(case_energies[i]) / (kT)))

      # including permutations term to account for probability
      Z_c_m_sums[e_i][t_i] = _accuracy(permutations[e_i] / _accuracy(cnt_case_energies)) * sum_Ediff
        
  return Z_c_m_sums