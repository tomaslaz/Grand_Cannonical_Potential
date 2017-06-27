"""
@author Tomas Lazauskas, 2016-2017

Input/output module.
"""

import numpy as np
import sys

import Constants
from source.Utilities import log

def check_file(file_path):
  """
  Checks if file exists.
  
  """
  
  return os.path.isfile(file_path)

def read_in_data(data_file, _accuracy, options):
  """
  Reads in all the necessary data from the input file.
  
  """
    
  names = None
  permutations = None
  chem_pot_multi = None
  data = None
  
  # reading in the energies
  log(__name__, "Reading in data: %s" % (data_file), options.verbose, indent=1)
  success, error, data = read_data(data_file, _accuracy)
  
  # reading in the permutations
  if success:
    success, error, permutations = read_line_values_as_array(data_file, _accuracy, line_no=2)
  
  # reading in chemical potential multipliers
  if success:
    success, error, chem_pot_multi = read_line_values_as_array(data_file, _accuracy, line_no=3)
      
  return success, error, names, permutations, chem_pot_multi, data

def read_line_values_as_array(file_path, _accuracy, line_no):
  """
  Reads in chemical potential multipliers
  
  """
  
  success = True
  error = ""
  values = None
  
  if not check_file:
    success = False
    error = "File: %s cannot be found." % (file_path)
  
  else:
    
    try:
      f = open(file_path)
      
    except:
      success = False
      error = "File: %s cannot be opened." % (file_path)
      return success, error, values
    
    line_cnt = 0
    for line in f:
      if line_cnt == (line_no-1):
        
        array = line.split(",")
        array_cnt = len(array)
              
        values = np.empty([array_cnt], dtype=_accuracy)
        
        for i in range(array_cnt):
          values[i] = _accuracy(array[i].strip())

      elif line_cnt > 1:
        break
      
      line_cnt += 1
      
    f.close()
  
  return success, error, values

def read_data(file_path, _accuracy):
  """
  Reads in the data
  
  """
  
  success = True
  error = ""
  energies_data = None
  
  if not check_file:
    success = False
    error = "File: %s cannot be found." % (file_path)
  
  else:
    energies_data = np.loadtxt(file_path, skiprows=3, delimiter=',', unpack=True, dtype=_accuracy)
    
  return success, error, energies_data

def read_permutations(file_path, _accuracy):
  """
  Reads in the permutations
  
  """
  
  success = True
  error = ""
  permutations = None
  
  if not check_file:
    success = False
    error = "File: %s cannot be found." % (file_path)
  
  else:
    
    try:
      f = open(file_path)
      
    except:
      success = False
      error = "File: %s cannot be opened." % (file_path)
      return success, error, permutations
    
    line_cnt = 0
    for line in f:
      if line_cnt == 1:
        
        array = line.split(",")
        array_cnt = len(array)
              
        permutations = np.empty([array_cnt], dtype=_accuracy)
        
        for i in range(array_cnt):
          permutations[i] = _accuracy(array[i].strip())

      elif line_cnt > 1:
        break
      
      line_cnt += 1
      
    f.close()
    
  return success, error, permutations

def write_Wm(temperatures, chem_pot_range, chem_pot_multi, Wm_array):
  """
  Writes the Wm array into a csv file
  
  """
  
  success = True
  error = ""
  
  try:
    fout = open(Constants.wm_results_filename, "w")
  except:
    success = False
    error = __name__ + ": Cannot open: " + filePath
    
    return success, error
  
  temp_len = len(temperatures)
  chem_pot_len = len(chem_pot_range)
  chem_pot_multi_len = len(chem_pot_multi)
  
  # first column is dedicated for the value of m
  fout.write(",")
  
  # writing mu values as header
  for mu_idx in range(chem_pot_len-1):
    fout.write("%f," % (chem_pot_range[mu_idx]))
    
  fout.write("%f\n" % (chem_pot_range[chem_pot_len-1]))
  
  # for each temperature:
  for t_i in range(temp_len):
    
    # writing temperature
    temperature = temperatures[t_i]
    fout.write("%f\n" % (temperature))
    
    # writing probabilities
    for m_index in range(chem_pot_multi_len):
      
      # writing the value of m
      fout.write("%f" % (chem_pot_multi[m_index]))
      
      for mu_idx in range(chem_pot_len):
        fout.write(",%e" % (Wm_array[t_i, mu_idx, m_index]))
      
      fout.write("\n")
          
  fout.close()
    
  return success, error