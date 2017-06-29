#!/usr/bin/python
"""
@author Tomas Lazauskas, 2016-2017

"""

import numpy as np
from optparse import OptionParser

import source.IO as IO
import source.GrandPotential as GrandPotential
from source.Utilities import log

_accuracy = np.float128

def cmd_line_args():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog input_data.csv"
  
  parser = OptionParser(usage=usage)

  parser.disable_interspersed_args()
  
  parser.add_option("-o", "--output", dest="outputData", action="store_true", default=False, 
    help="A flag to write the data into files")
  parser.add_option("-p", "--permutations", dest="permCalc", action="store_true", default=False, 
    help="Permutations will be calculated instead of read from the data file. User must modify: DistributionAnalysis.calc_permutation")
  parser.add_option('-t', dest="temps", default="1", help="List of temperatures separated by a comma (default: T=0)")
  parser.add_option('-u', dest="urange", default=None, 
                    help="Chemical potential range: (default: None, syntax: 0,0.1,0.1)")
  parser.add_option("-v", dest="verbose", default=1, type="int", help="Verbose: 0 - off, 1 - on.")
  
  (options, args) = parser.parse_args()

  if (len(args) != 1):
    parser.error("incorrect number of arguments")

  return options, args

def main(options, args):
  """
  The main routine.
  
  """
  
  data_file = args[0]
  
  # reading in the data
  success, error, names, permutations, chem_pot_multi, data = IO.read_in_data(data_file, _accuracy, options)
  
  if success:
    # grand canonical analysis
    log(__name__, "Performing the grand canonical analysis", options.verbose, indent=1)
    success, error = GrandPotential.perform_grand_canonical_analysis(names, permutations, chem_pot_multi, data, _accuracy, options)
  
  return success, error

if __name__ == "__main__":
  
  # command line arguments and options
  options, args = cmd_line_args()
  
  log(__name__, "Grand Canonical Potential", options.verbose)
  
  # the main routine
  success, error = main(options, args)
  
  if success:
    log(__name__, "Finished.", options.verbose)
    
  else:
    log(__name__, "ERROR: %s" % (error), 1)
  