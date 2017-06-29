"""
@author Tomas Lazauskas, 2016-2017

Utilities module.
"""

import datetime
import numpy as np

def frange(x, y, jump):
  frange_arr = []
  
  while x < y:
    frange_arr.append(x)
    x += jump
  
  return frange_arr

def factorial_division(up_fac, down_fac, _accuracy):
  """
  Calculates a division of two factorials
  
  """
  
  up_fac_int = np.int32(up_fac)
  down_fac_int = np.int32(down_fac)
  
  result = _accuracy(1)
  
  if up_fac_int > down_fac_int:
    for i in range(down_fac_int+1, up_fac_int+1):
      result *= _accuracy(i)
      
  elif up_fac_int < down_fac_int:
    for i in range(up_fac_int+1, down_fac_int+1):
      result *= _accuracy(i)
    
    result = _accuracy(1) / result
    
  return result

def get_chem_pot_range(chem_pot_string):
  """
  Get the list of chemical potential values for the analysis
  
  """
  
  success = True
  error = ""
  chem_pot_range = None
  
  chem_pot_range_arr = chem_pot_string.split(",")
  
  if len(chem_pot_range_arr) != 3:
    success = False
  
  chem_pot_range = np.arange(float(chem_pot_range_arr[0]), float(chem_pot_range_arr[1]), float(chem_pot_range_arr[2]))
    
  if len(chem_pot_range) < 1:
    success = False

  if not success:
    error = "Incorrect chemical potential range: %s" % (chem_pot_string)
    return success, error, None
  
  return success, error, chem_pot_range

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

def get_temperature_colour(temperature):
  """
  Returns colour name with respect to the temperature
  
  """
  
  colour = 'black'
  
  if temperature <= 100:
    colour = 'b'
      
  elif temperature <= 300:
    colour = 'g'
    
  elif temperature <= 500:
    colour = 'r'
    
  elif temperature <= 1000:
    colour = 'c'
  
  elif temperature <= 2000:
    colour = 'y'

  elif temperature <= 3000:
    colour = 'm'

  elif temperature <= 4000:
    colour = 'darkblue'

  elif temperature <= 5000:
    colour = 'sienna'

  elif temperature <= 6000:
    colour = 'indigo'
  
  elif temperature <= 7000:
    colour = 'orange'
  
  elif temperature <= 8000:
    colour = 'grey'
    
  elif temperature <= 9000:
    colour = 'brown'
    
  elif temperature <= 10000:
    colour = 'black'
  
  return colour

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