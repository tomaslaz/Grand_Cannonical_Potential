"""
@author Tomas Lazauskas, 2017

A module to analyse distribution functions of different stoichiometries

"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import Constants
import Utilities
from Utilities import log

def calc_permutation(m, mm, _accuracy):
  """
  Evaluates the permutation expression
  
  """
  
  cat_sites = 108
  an_sites = 216
  
  up_fac = cat_sites - 2*m
  down_fac = cat_sites - 2*mm
  
  value = Utilities.factorial_devision(up_fac, down_fac, _accuracy)
  
  value *= Utilities.factorial_devision(2*m, 2*mm, _accuracy)
  
  up_fac = an_sites - m
  down_fac = an_sites - mm
  
  value *= Utilities.factorial_devision(up_fac, down_fac, _accuracy)
  
  value *= Utilities.factorial_devision(m, mm, _accuracy)
  
  return value

def distribution_analysis(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, 
                        experiment_cnts, _accuracy, options):
  """
  Performs the distribution analysis: evaluates Wm and plots it against m and mu.
  
  """
  
  success = True
  error = ""
  
  # Preparing Wm probabilities    
  Wm_array = prepare_Wm(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, 
                        experiment_cnts, _accuracy, options)
  
  # Plotting the Wm probabilities
  dist_func_analysis_temp_mu_3D_plots(temperatures, chem_pot_range, chem_pot_multi, Wm_array, _accuracy, options)
  
  return success, error

def dist_func_analysis_temp_mu_3D_plots(temperatures, chem_pot_range, chem_pot_multi, Wm_array, _accuracy, options):
  """
  Plots distribution functions in terms of m and mu with respect to the temperature
  
  """
  
  #print "chem_pot_range: ", chem_pot_range
#   print ": ", 
#   print ": ", 
  
  temp_len = len(temperatures)
  chem_pot_len = len(chem_pot_range)
  
  # for each temperature:
  for t_i in range(temp_len):
    
    temperature = temperatures[t_i]
    
    file_name = 'wm_%.2d.png' % (temperature)
    
    log(__name__, "Plotting the distribution functions for %.2d K (%s)" % (temperature, file_name), options.verbose, indent=3)
    
    fig = plt.figure(figsize=(Constants.fig_size_x, Constants.fig_size_y))
    ax = fig.add_subplot(111, projection='3d')
    
    # plotting for each mu
    for mu_i in range(chem_pot_len):
      mu = chem_pot_range[mu_i]
      
      x = chem_pot_multi
      y = Wm_array[t_i, mu_i, :]

      colour = Utilities.get_temperature_colour(temperature)
      
      if mu_i == 0:
        ax.plot(xs=x, ys=y, zs=mu, zdir='y', color=colour, linewidth=1.0, label="%d" % (temperature))
      else:
        ax.plot(xs=x, ys=y, zs=mu, zdir='y', color=colour, linewidth=1.0)
      
    ax.set_xlabel(r'$m$', fontsize=Constants.fig_label_fontsize)
         
    ax.set_zlabel(r'$w_{m}$', fontsize=Constants.fig_label_fontsize)
    ax.set_zlim3d(0, 1)
      
    ax.set_ylabel(r'$\mu (eV)$', fontsize=Constants.fig_label_fontsize)
     
    ax.set_ylim3d(np.min(chem_pot_range), np.max(chem_pot_range))
     
    plt.legend(loc=2, title="Temperature (K)", ncol=2, fontsize=Constants.fig_legend_fontsize, borderaxespad=0.)
     
    fig.savefig(file_name, dpi=Constants.fig_dpi, bbox_inches='tight')
     
    # clearing the figure settings   
    plt.clf()
    plt.close()
    
def prepare_Wm(chem_pot_multi, temperatures, chem_pot_range, min_energies, delta_E_sums, experiment_cnts, 
               _accuracy, options):
  """
  Evaluates Wm with respect to temperature and chemical potential
  
  """
  
  log(__name__, "Preparing Wm", options.verbose, indent=3)
  
  temp_len = len(temperatures)
  chem_pot_len = len(chem_pot_range)
  chem_pot_multi_len = len(chem_pot_multi)
  
  Wm_array = np.empty([temp_len, chem_pot_len, chem_pot_multi_len], _accuracy)
  
  for t_i in range(temp_len):
    kT = np.longdouble(Constants.kB * temperatures[t_i])
  
    for mu_i in range(chem_pot_len):    
      mu = chem_pot_range[mu_i]
      
      Wm_sum = _accuracy(0.0)
      for m_idx in range(chem_pot_multi_len):
        m_value = chem_pot_multi[m_idx]
        Emin_m = min_energies[m_idx]
        
        # Estimating the top part of Wm
        sum_top = delta_E_sums[m_idx][t_i]
        
        # Estimating the bottom part of Wm
        sum_bottom = _accuracy(0.0)
        
        for mm_index in range(chem_pot_multi_len):
          mm_value = chem_pot_multi[mm_index]
          Emin_mm = min_energies[mm_index]
          
          # exponential term
          exp_expr_pow = _accuracy(((Emin_m - Emin_mm) + mu*(m_value - mm_value))/kT)
          exp_expr = _accuracy(np.exp(exp_expr_pow))
          
          # Nm/Nmm
          attempts_cnt = experiment_cnts[m_idx] / experiment_cnts[mm_index] 
          
          # Pmm/Pm
          permutation = _accuracy(calc_permutation(m_value, mm_value, _accuracy))

          sum_bottom += exp_expr * attempts_cnt * permutation * delta_E_sums[mm_index][t_i]
        
        Wm_value = sum_top / sum_bottom
        
        Wm_sum += Wm_value
        
        Wm_array[t_i, mu_i, m_idx] = Wm_value
      
      # Normalising
      Wm_array[t_i, mu_i, :] /= Wm_sum
      
  return Wm_array