# Grand_Cannonical_Potential

The work is based on the methodology developed in

Usage: GrandPotential.py input_data.csv

Options:
  -h, --help          show this help message and exit
  -o, --output        A flag to write the data into files
  -p, --permutations  Permutations will be calculated instead of read from the
                      data file. User must modify:
                      DistributionAnalysis.calc_permutation
  -t TEMPS            List of temperatures separated by a comma (default: T=0)
  -u URANGE           Chemical potential range: (default: None, syntax:
                      0,0.1,0.1)
  -v VERBOSE          Verbose: 0 - off, 1 - on.

input_data file must be in csv format and the energies should be given in columns where: 
1) the first line is the name of the simulation, 
2) the second line is the number of total possible configurations 
3) the third line is an integer value to be used with a chemical potential. If the chemical potential is not used, the third line should be left blank.

Important: all columns in the input_data.csv file should have the same number of lines. The missing values should be replaced with 0 which will not be included in the calculations.
