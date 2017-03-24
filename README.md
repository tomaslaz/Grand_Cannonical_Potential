# Grand_Cannonical_Potential

The work is based on the methodology developed in

Usage: GrandPotential.py input_data.csv

Options:
  -h, --help  show this help message and exit
  -t TEMPS    List of temperatures separated by a comma (default t=0)
  -v VERBOSE  Verbose: 0 - off, 1 - on (default v=1).

input_data file must be in csv format and the energies should be given in columns where: 1) the first line is the name of the simulation, 2) the second line is the number of total possible configurations 3) the third line is an integer value to be used with a chemical potential. If the chemical potential is not used, the third line should be left blank.

Important: all columns in the input_data.csv file should have the same number of lines. The missing values should be replaced with 0 which will not be included in the calculations.
