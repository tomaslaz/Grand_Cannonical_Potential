"""
@author Tomas Lazauskas, 2017

A module for plotting the results

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import Constants
import Utilities
from Utilities import log

def _setup_temperature_legend(ax):
  """
  Sets up a legend for temperatures
  
  """
  
  lg = ax.legend(loc=3, bbox_to_anchor=(0.,1.02,1.,0.102),
                 title="Temperature (K)", ncol=5, fontsize=Constants.fig_legend_fontsize, mode="expand", borderaxespad=0.)
  
  lg.get_frame().set_alpha(0.5)
  lg.get_title().set_fontsize('%d' % (Constants.fig_label_fontsize))
  
def avg_values(temperatures, chem_pot_range, avg_array, prop_name, _accuracy, options):
  """
  Plots average values with respect with mu and temperature
  
  """
  
  fig = plt.figure(figsize=(Constants.fig_size_x, Constants.fig_size_y))
  ax = fig.add_subplot(1,1,1)
  
  temp_len = len(temperatures)
  
  # for each temperature:
  for t_index in range(temp_len):
    temperature = temperatures[t_index]
    
    colour = Utilities.get_temperature_colour(temperature)
    
    ax.plot(chem_pot_range, avg_array[t_index, :], color=colour, linewidth=1.5, label="%.2d" % (temperature))
  
  plt.grid()
  
  ax.set_xlabel(r'$\mu$', fontsize = Constants.fig_label_fontsize)
  ax.set_ylabel(r'$<%s>$' % (prop_name), fontsize=Constants.fig_label_fontsize)
  
  # place the legend
  _setup_temperature_legend(ax)
  
  fig.savefig(Constants.avg_plot_filename % (prop_name), dpi=Constants.fig_dpi, bbox_inches='tight')
  
  # clearing the figure settings   
  plt.clf()
  plt.close()

def c_g_omega_over_mu(temperatures, chem_pot_range, omega_arr, _accuracy, options):
  """
  Plots omega with respect with mu and temperature
  
  """
  
  fig = plt.figure(figsize=(Constants.fig_size_x, Constants.fig_size_y))
  ax = fig.add_subplot(1,1,1)
  
  temp_len = len(temperatures)
  
  # for each temperature:
  for t_index in range(temp_len):
    temperature = temperatures[t_index]
  
    colour = Utilities.get_temperature_colour(temperature)
    
    ax.plot(chem_pot_range, omega_arr[t_index, :], color=colour, linewidth=1.5, label="%.2d" % (temperature))
    
  plt.grid()
  
  ax.set_xlabel(r'$\mu$', fontsize = Constants.fig_label_fontsize)
  ax.set_ylabel(r'$\Omega$', fontsize=Constants.fig_label_fontsize)
  
   # place the legend
  _setup_temperature_legend(ax)
  
  fig.savefig(Constants.omega_mu_plot_filename, dpi=Constants.fig_dpi, bbox_inches='tight')
  
  # clearing the figure settings   
  plt.clf()
  plt.close()

def c_omega(chem_pot_multi, names, temperatures, omega_arr, _accuracy, options):
  """
  """
  
  fig = plt.figure(figsize=(Constants.fig_size_x, Constants.fig_size_y))
  ax = fig.add_subplot(1,1,1)
  
  temp_len = len(temperatures)
  
  # for each temperature:
  for t_index in range(temp_len):
    temperature = temperatures[t_index]
  
    colour = Utilities.get_temperature_colour(temperature)
  
    ax.plot(chem_pot_multi, omega_arr[t_index, :], color=colour, linewidth=1.5, label="%.2d" % (temperature))
  
  plt.grid()
  
  ax.set_xlabel(r'$m$', fontsize = Constants.fig_label_fontsize)
  ax.set_ylabel(r'$\gamma^{c}_{m}$', fontsize=Constants.fig_label_fontsize)
  
  # x axis ticklabels
  plt.xticks(chem_pot_multi, names, rotation=30, fontsize = Constants.fig_ticklabel_fontsize)
  
  # place the legend
  _setup_temperature_legend(ax)
  
  fig.savefig(Constants.omega_plot_filename, dpi=Constants.fig_dpi, bbox_inches='tight')
  
  # clearing the figure settings   
  plt.clf()
  plt.close()

def wm_contour(temperatures, chem_pot_range, chem_pot_multi, Wm_array, _accuracy, options):
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
  
    contour_filled = plt.contourf(X, Y, Z, levels, cmap=cmap, interpolation='bilinear', vmax=1.1, vmin=0.0, alpha=0.95)
    
    cbar = plt.colorbar(contour_filled, ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    cbar.ax.tick_params(labelsize=16) 
    cbar.ax.set_ylabel(r'$w_{m}$', fontsize=Constants.fig_label_fontsize)
      
    for tick in cbar.ax.yaxis.get_major_ticks():
      tick.label.set_fontsize(16) 
      
    plt.xlabel(r'$m$', fontsize = Constants.fig_label_fontsize)
    plt.ylabel(r'$\mu$', fontsize = Constants.fig_label_fontsize)
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