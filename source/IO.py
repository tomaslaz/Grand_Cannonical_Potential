"""
@author Tomas Lazauskas, 2016-2017

Input/output module.
"""

import numpy as np
import sys

def check_file(file_path):
  """
  Checks if file exists.
  
  """
  
  return os.path.isfile(file_path)

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
    energies_data = np.loadtxt(file_path, skiprows=2, delimiter=',', unpack=True, dtype=_accuracy)
    
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