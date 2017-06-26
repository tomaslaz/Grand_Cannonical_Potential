"""
@author Tomas Lazauskas, 2017

A module to analyse distribution functions of different stoichiometries

"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import Constants
import sys
import Utilities
from Utilities import log
import warnings

# TODO: needs to made universal
def calc_permutation(m, mm, _accuracy):
  """
  Evaluates the permutation expression
  
  """
  
  cat_sites = 108
  an_sites = 216
  
  up_fac = cat_sites - 2*m
  down_fac = cat_sites - 2*mm
  
  value = Utilities.factorial_division(up_fac, down_fac, _accuracy)
  
  value *= Utilities.factorial_division(2*m, 2*mm, _accuracy)
  
  up_fac = an_sites - m
  down_fac = an_sites - mm
  
  value *= Utilities.factorial_division(up_fac, down_fac, _accuracy)
  
  value *= Utilities.factorial_division(m, mm, _accuracy)
  
  return value

def distribution_analysis(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, 
                        experiment_cnts, permutations, _accuracy, options):
  """
  Performs the distribution analysis: evaluates Wm and plots it against m and mu.
  
  """
  
  success = True
  error = ""
  
  # Preparing Wm probabilities    
  Wm_array = prepare_Wm(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, 
                        experiment_cnts, permutations, _accuracy, options)
  
  # Plotting the Wm probabilities 3D plots
  dist_func_analysis_temp_mu_m_contour(temperatures, chem_pot_range, chem_pot_multi, Wm_array, _accuracy, options)
  
  return success, error

def dist_func_analysis_temp_mu_m_contour(temperatures, chem_pot_range, chem_pot_multi, Wm_array, _accuracy, options):
  """
  Plots distribution functions in terms of m and mu with respect to the temperature
  
  """
  
  cmap = plt.cm.OrRd
  
  levels = np.arange(0.0, 1.05, 0.05)
  
  temp_len = len(temperatures)
  chem_pot_len = len(chem_pot_range)
  
  X, Y = np.meshgrid(chem_pot_multi, chem_pot_range)
  
  # for each temperature:
  for t_i in range(temp_len):
    
    temperature = temperatures[t_i]
    
    file_name = 'wm_%.2d.png' % (temperature)
    
    log(__name__, "Plotting distribution functions @ %.2d K (%s)" % (temperature, file_name), options.verbose, indent=3)
    
    fig = plt.figure(figsize=(Constants.fig_size_x, Constants.fig_size_y))

    Z = Wm_array[t_i, :, :]    
  
    contour_filled = plt.contourf(Y, X, Z, levels, cmap=cmap, interpolation='bilinear', vmax=1.1, vmin=0.0)
    
    cbar = plt.colorbar(contour_filled, ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    cbar.ax.tick_params(labelsize=16) 
    cbar.ax.set_ylabel(r'$w_{m}$', fontsize=Constants.fig_label_fontsize)
      
    for tick in cbar.ax.yaxis.get_major_ticks():
      tick.label.set_fontsize(16) 
      
      
    plt.xlabel(r'$\mu$', fontsize = Constants.fig_label_fontsize)
    
    plt.ylabel(r'$<m>$', fontsize = Constants.fig_label_fontsize)
      
    plt.title("T = %d K" % (temperature), fontsize = Constants.fig_legend_fontsize)
    
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
      tick.label.set_fontsize(16) 
      
    for tick in ax.yaxis.get_major_ticks():
      tick.label.set_fontsize(16)
    
    fig.savefig(file_name+".png", dpi=Constants.fig_dpi, bbox_inches='tight')
         
    # clearing the figure settings   
    plt.clf()
    plt.close()
    
def prepare_Wm(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, experiment_cnts, 
               permutations, _accuracy, options):
  """
  Evaluates Wm with respect to temperature and chemical potential
  
  """
  
  log(__name__, "Preparing Wm", options.verbose, indent=3)
  
  temp_len = len(temperatures)
  chem_pot_len = len(chem_pot_range)
  chem_pot_multi_len = len(chem_pot_multi)
  
  Wm_array = np.zeros([temp_len, chem_pot_len, chem_pot_multi_len], _accuracy)
  
  for t_i in range(temp_len):
    kT = np.longdouble(Constants.kB * temperatures[t_i])
  
    for mu_i in range(chem_pot_len):    
      mu = chem_pot_range[mu_i]
      
      breakLoops = False
      
      Wm_sum = _accuracy(0.0)
      for m_index in range(chem_pot_multi_len):
        
        m_value = chem_pot_multi[m_index]
        Emin_m = min_energies[m_index]
        
        # Estimating the top part of Wm
        sum_top = delta_E_sums[m_index][t_i]
        
        # Estimating the bottom part of Wm
        sum_bottom = _accuracy(0.0)
        
        # Calculating exponential terms powers
        exp_powers = np.empty(chem_pot_multi_len, _accuracy)
        
        for mm_index in range(chem_pot_multi_len):
          mm_value = chem_pot_multi[mm_index]
          Emin_mm = min_energies[mm_index]
          
          # the power of the exponential term
          exp_powers[mm_index] = _accuracy(((Emin_m - Emin_mm) + mu*(m_value - mm_value))/kT)
        
        for mm_index in range(chem_pot_multi_len):
          if mm_index != m_index:
            
            mm_value = chem_pot_multi[mm_index]
            Emin_mm = min_energies[mm_index]
            
            # exponential term
            exp_expr_pow = exp_powers[mm_index]
                        
            if (exp_expr_pow > Constants.BIGEXPO):
              # the value of the exponential power is too damn high
              breakLoops = True              
              break
            
            elif exp_expr_pow < -Constants.BIGEXPO:
              exp_expr = 0.0
            else:
               exp_expr = _accuracy(np.exp(exp_expr_pow))
            
            # Nm/Nmm
            attempts_cnt = experiment_cnts[m_index] / experiment_cnts[mm_index] 
            
            # Pmm/Pm
            if options.permCalc:
              permutation = _accuracy(calc_permutation(m_value, mm_value, _accuracy))
            else:
              permutation = _accuracy(_accuracy(permutations[mm_index])/_accuracy(permutations[m_index]))         
            
            sum_bottom += exp_expr * attempts_cnt * permutation * delta_E_sums[mm_index][t_i]
          
        if not breakLoops:
          
          if sum_bottom == 0.0:
            Wm_value = 1.1
          else:
            Wm_value = sum_top / sum_bottom
            
          Wm_sum += Wm_value
        else:
          Wm_value = 0

        Wm_array[t_i, mu_i, m_index] = Wm_value
       
      # Normalising
      if Wm_sum != 0.0:
        Wm_array[t_i, mu_i, :] /= Wm_sum
  
  return Wm_array