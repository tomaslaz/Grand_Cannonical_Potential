# Grand_Cannonical_Potential

Usage: GrandPotential.py input_data.csv

Options:
  -h, --help  show this help message and exit
  -t TEMPS    List of temperatures separated by a comma (default t=0)
  -v VERBOSE  Verbose: 0 - off, 1 - on.

input_data file must be in csv format and the energies should be given in columns where the first line is the name of the simulation and the second line is the number of total possible configurations. 

Important: all columns in the input_data.csv file should have the same number of lines. The missing values should be replaced with 0 which will not be included in the calculations.