#!/usr/bin/env python

#Plotting parameters
#Tural Aksel

import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as py
import matplotlib.cm as cm
import numpy as np
import math

params = {'axes.labelsize':40,
          'axes.titlesize':40,
          'font.size': 40,
          'legend.fontsize': 40,
          'legend.labelspacing':0.1,
          'legend.handletextpad':0.2,
          'legend.borderaxespad':1.0,
          'legend.numpoints':1,
          'legend.handlelength':0.5,
          'legend.frameon'  : True,
          'legend.fancybox' : True,
          'legend.borderpad': 0.2,
          'lines.markersize': 25,
          'lines.markeredgewidth': 3,
          'axes.linewidth'  : 1.0,
          'lines.linewidth' : 10,
          'mathtext.fontset':'cm',
          'mathtext.default' :'regular',
          'mathtext.fallback_to_cm' : True,
          'axes.formatter.limits' : (-4, 4),
          'figure.subplot.top'    : 0.95,
          'figure.subplot.bottom' : 0.135
          }

py.rcParams.update(params)

font = {'family' : 'sans-serif',
        'sans-serif':'Arial',
        'weight' : 'normal',
        }

py.rc('font', **font)  # pass in the font dict as kwargs


yTick = {'major.pad': 10,'major.size': 20,'minor.size': 10,'labelsize':40}
xTick = {'major.pad': 10,'major.size': 20,'minor.size': 10,'labelsize':40}

py.rc('xtick', **xTick)  # pass in the font dict as kwargs
py.rc('ytick', **yTick)  # pass in the font dict as kwargs

#Convert between pixels -> inches / find the best ratio image using golden ratio
def get_figsize(fig_width_pt):
    inches_per_pt = 1.0/72.0                # Convert pt to inch
    golden_mean = (math.sqrt(5)-1.0)/2.0    # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    fig_size =  [fig_width,fig_height]      # exact figsize
    return fig_size

def make_N_colors(cmap_name, N): 
    cmap = cm.get_cmap(cmap_name, N) 
    return cmap(np.arange(N)) 
