#!/usr/bin/python
"""
@author Tomas Lazauskas, 2016-2017

"""

#import copy
#import math
#import numpy as np
from optparse import OptionParser
#import sys

from source.Utilities import log

def cmd_line_args():
  """
  Handles command line arguments and options.
  
  """
  
  usage = "usage: %prog"
  
  parser = OptionParser(usage=usage)

  parser.disable_interspersed_args()
  
  parser.add_option('-t', dest="temps", default=None, help="List of temperatures separated by a comma (default t=0)")
  parser.add_option("-v", dest="verbose", default=0, type="int", help="Verbose: 0 - off, 1 - on.")
  
  (options, args) = parser.parse_args()

  if (len(args) != 0):
    parser.error("incorrect number of arguments")

  return options, args

def main(options, args):
  """
  
  """
  
  print ""

if __name__ == "__main__":
  
  # command line arguments and options
  options, args = cmd_line_args()
  
  log(__name__, "Grand Cannonical Potential", options.verbose)
  
  main(options, args)
  
  log(__name__, "Finished.", options.verbose)