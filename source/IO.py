"""
@author Tomas Lazauskas, 2016-2017

Input/output module.
"""

import numpy as np

def check_file(file_path):
  """
  Checks if file exists.
  
  """
  
  return os.path.isfile(file_path)

def read_data(file_path):
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
    energies_data = np.loadtxt(file_path, skiprows=1, delimiter=',', unpack=True)
  
  return success, error, energies_data