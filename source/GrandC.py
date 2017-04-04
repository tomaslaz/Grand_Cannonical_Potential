"""
@author Tomas Lazauskas, 2016-2017

Grand canonical ensemble module.
"""

import copy
import math
import numpy as np
import sys

import Constants
import DistributionAnalysis
import Utilities
from Utilities import log
  
def perform_grand_canonical_analysis(names, permutations, chem_pot_multi, data, _accuracy, options):
  """
  The main method of the grand canonical analysis
  
  """
  
  chem_pot_analysis = False
  
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
    
    if not success:
      return success, error
    
    chem_pot_analysis = True
      
  log(__name__, "Analysis will be performed at (K): %s" % (options.temps), options.verbose, indent=2)
  
  # first of all lets prepare the data for the calculations.
  log(__name__, "Preparing the data", options.verbose, indent=2)
  
  energies, min_energies, shifted_energies, experiment_cnts = prepare_energies(data, _accuracy, options)
  
  delta_E_sums = prepare_delta_E_sums(shifted_energies, temperatures, _accuracy, options)
  
  # Perform chemical potential analysis
  if chem_pot_analysis:
    # performs distribution analysis    
    success, error = DistributionAnalysis.distribution_analysis(chem_pot_multi, temperatures, chem_pot_range, 
                                           min_energies, delta_E_sums, experiment_cnts, _accuracy, options)
      
  if not success:
    return success, error
  
  return success, error

def prepare_delta_E_sums(shifted_energies, temperatures, _accuracy, options):
  """
  Prepares the sum of the second part of the Z_c_m equation
  
  """
    
  energies_cnt = len(shifted_energies)
  temp_cnt = len(temperatures)
  
  delta_E_sums = np.empty([energies_cnt, temp_cnt], dtype=_accuracy)
    
  log(__name__, "Preparing delta_E_sums to speed up the calculations", options.verbose, indent=3)
  
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

      delta_E_sums[e_i][t_i] = sum_Ediff
    
  return delta_E_sums

def prepare_energies(input_data_array, _accuracy, options):
  """
  Prepares energy arrays for the calculations
  
  """
  
  log(__name__, "Preparing the energies", options.verbose, indent=3)
  
  cases_cnt = input_data_array.shape[0]
  
  min_energies = np.empty(cases_cnt, dtype=_accuracy)
  energies = []
  shifted_energies = []
  experiment_cnts = np.empty(cases_cnt, dtype=_accuracy)
  
  for i in range(cases_cnt):
    # removes all zero values from the beginning and end of the array
    case_energies = np.trim_zeros(copy.deepcopy(input_data_array[i]))
    energies.append(case_energies)
    
    # it might be the first value but it might not be, thus we getting the minimum
    min_energies[i] = _accuracy(np.min(case_energies))
    
    diff_energy = case_energies - min_energies[i]
    shifted_energies.append(diff_energy)
    
    # saving the number of experiments of the stoichiometry
    experiment_cnts[i] = len(case_energies)
  
  return energies, min_energies, shifted_energies, experiment_cnts

  