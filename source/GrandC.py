"""
@author Tomas Lazauskas, 2016-2017

Grand canonical ensemble module.
"""

import copy
import numpy as np

from Utilities import log

def perform_grand_canonical_analysis(input_data_array, temperatures, _accuracy, options):
  """
  The main method of the grand canonical analysis
  
  """
  
  success = True
  error = ""
  
  energy_list = []
  name_list = []
  total_combi_list = []
  
  # first of all lets prepare the data for the calculations.
  log(__name__, "Preparing the data", options.verbose, indent=2)
  
  energies, min_energies, shifted_energies = prepare_energies(input_data_array, _accuracy, options)
  
  
  
  return success, error

def prepare_energies(input_data_array, _accuracy, options):
  """
  Prepares energy arrays for the calculations
  
  """
  
  log(__name__, "Preparing the energies", options.verbose, indent=3)
  
  cases_cnt = input_data_array.shape[0]
  
  min_energies = np.zeros(cases_cnt, _accuracy)
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