"""
@author Tomas Lazauskas, 2016-2017

Utilities module.
"""

import datetime

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