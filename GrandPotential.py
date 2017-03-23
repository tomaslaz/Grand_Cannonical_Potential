#!/usr/bin/python
"""
@author Tomas Lazauskas, 2016-2017

"""

#import copy
#import math
#import numpy as np
from optparse import OptionParser
#import sys

import source.IO as IO
import source.GrandC as GrandC
from source.Utilities import log

def cmd_line_args():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog input_data.csv"
  
  parser = OptionParser(usage=usage)

  parser.disable_interspersed_args()
  
  parser.add_option('-t', dest="temps", default=None, help="List of temperatures separated by a comma (default t=0)")
  parser.add_option("-v", dest="verbose", default=0, type="int", help="Verbose: 0 - off, 1 - on.")
  
  (options, args) = parser.parse_args()

  if (len(args) != 1):
    parser.error("incorrect number of arguments")

  return options, args

def main(options, args):
  """
  The main routine.
  
  """
  
  # reading in the data
  success, error, data = IO.read_data(args[0])
  
  if success:
    # grand cannonical analysis
    success, error = GrandC.perform_grand_canonical_analysis(data)
  
  return success, error
  
if __name__ == "__main__":
  
  # command line arguments and options
  options, args = cmd_line_args()
  
  log(__name__, "Grand Cannonical Potential", options.verbose)
  
  # the main routine
  main(options, args)
  
  log(__name__, "Finished.", options.verbose)
  