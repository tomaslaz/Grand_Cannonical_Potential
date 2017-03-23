"""
@author Tomas Lazauskas, 2016-2017

Utilities module.
"""

import datetime
import numpy as np

def get_list_of_temps(temp_string):
  """
  A function to process an argument string line to an array of temperatures
  
  """
  success = True
  error = ""
  temps = None
  
  arr_string = temp_string.split(",")
  arr_string_len = len(arr_string)
  
  temps = np.zeros(arr_string_len)
  
  for i in range(arr_string_len):
    try:
      temp = float(arr_string[i])
    except:
      success = False
      error = "Incorrect temperatures: %s" % (temp_string)
      break
    
    temps[i] = temp
  
  return success, error, temps

def log(caller, message, verbose, indent=0):
  """
  Log output to screen
    
  indent:
      how much to indent the output by (defaults to zero)
  
  """
  
  if verbose:
    now = datetime.datetime.now().strftime("%d/%m/%y, %H:%M:%S: %f")
    ind = ""
    for _ in xrange(indent):
        ind += "  "
        
    print "[%s]: %s%s >> %s" % (now, ind, caller, message)