# -*- coding: utf-8 -*-

"""
This script compiles the C libraries.

You should ensure this script compiles automatically on every 
machine you intend to run on.

"""
import os
import sys
import glob

# guess correct settings
# Mac?
if os.uname()[0] == "Darwin":
  # use gcc
  CC = "gcc "
  CCFLAGS = "-c -O3 -fPIC "
  CCLINKFLAGS = "-bundle -flat_namespace -undefined suppress "  
  INCLUDES = ""

# Linux box?
elif os.uname()[0] == "Linux":

  CC = "gcc "
  CCFLAGS = "-c -O3 -fPIC "
  CCLINKFLAGS = "-shared "  
  INCLUDES = ""

################################################################################
# NOTHING BELOW THIS LINE SHOULD NEED TO BE CHANGED
################################################################################

if os.path.dirname(__file__):
    os.chdir(os.path.dirname(__file__))


OWD = os.getcwd()
os.chdir("source")
os.chdir("clibs")

# MAKE SURE NO WRAP C LIBS LEFT OVER
wrap_c_files = glob.glob("*_wrap.c")
if len(wrap_c_files):
  for fn in wrap_c_files:
    os.unlink(fn)

CLIBSRC = [fn[:-2] + "c" for fn in glob.glob("*[!__init__].py") if not fn.endswith("_utils.py")]

DEPS = {}

SOBJ = {}
OBJ = {}
for foo in CLIBSRC:
  SOBJ[foo] = "_%s.so" % foo[:-2]
  OBJ[foo] = "%s.o" % foo[:-2]

EXTC = []
tmpc = glob.glob("*.c")
for foo in tmpc:
  match = False
  if foo in CLIBSRC:
    match = True
  
  if match:
    pass
  
  else:
    EXTC.append(foo)

def checkForExe(exe):
  """
  Check if executable can be located 
  
  """
  # check if exe programme located
  syspath = os.getenv("PATH", "")
  syspatharray = syspath.split(":")
  found = 0
  for syspath in syspatharray:
    if os.path.exists(os.path.join(syspath, exe)):
      found = 1
      break
  
  if found:
    exepath = exe
  
  else:
    exepath = 0
  
  return exepath

def runCommand(command, checkStatus=1):
  
  status = os.system(command)
  if checkStatus and status:
    print "COMMAND FAILED"
    sys.exit(25)

def main():
    
  os.chdir(OWD)
  
  print "cd source"
  os.chdir("source")
  print "cd clibs"
  os.chdir("clibs")
  print "-------------"
  
  # first compile non swig C libraries
  for source in EXTC:
    command = "%s %s %s %s" % (CC, CCFLAGS, INCLUDES, source)
    print command
    runCommand(command)

  print "-------------"
  
  # compile non swig c files
  for src in CLIBSRC:
    # compile src file
    command = "%s %s %s" % (CC, CCFLAGS, src)
    print command
    runCommand(command)
    
    print "-------------"
  
  # link non swig libraries
  for src in CLIBSRC:
    try:
      EXTO = DEPS[src]
        
    except KeyError:
      EXTO = ""
        
    command = "%s %s %s %s -o %s" % (CC, CCLINKFLAGS, OBJ[src], EXTO, SOBJ[src])
    print command
    runCommand(command)
    
    print "-------------"
    
def clean():
    
  os.chdir(OWD)
  
  print "cd source"
  os.chdir("source")
  
  print "cd clibs"
  os.chdir("clibs")
  
  # remove compiled files
  command = "rm -f *.pyc *.o *.so *_wrap.c "
  print command
  runCommand(command, checkStatus=0)
  print "-------------"
  
  os.chdir("..")
  
  # remove python compiled files
  command = "rm -f *.pyc"
  print command
  runCommand(command, checkStatus=0)
  print "-------------"

if __name__ == '__main__':
  if len(sys.argv) > 1:
    if sys.argv[1] == "clean":
      clean()
      sys.exit(0)
      
    elif sys.argv[1] == "all":
      pass
    
    else:
      print "Usage: setup.py [options]\n\nOptions are:\n    clean"
      sys.exit(8)
  
  clean()
  main()
