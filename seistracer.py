#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
 SeisTracer

 Author: Caio Ciardelli, Northwestern University, February 2023

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

-----------------------------------------------------------------------------------------------
"""

# Import necessary modules
import os
import sys
import copy
import ctypes as cp
import math
import struct
import numpy as np
import warnings
import matplotlib
import matplotlib.pyplot as plt

# Store default matplotlib backend
DEFAULT_BACKEND = matplotlib.get_backend ()

# Import utility modules
from glob import glob
from pathlib import Path
from scipy.interpolate import interp1d
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

# Load shared libraries
libcoefficients = cp.cdll.LoadLibrary (os.path.dirname (os.path.realpath (__file__)) + '/lib/libcoefficients.so')
libio           = cp.cdll.LoadLibrary (os.path.dirname (os.path.realpath (__file__)) + '/lib/libio.so')
libtracer       = cp.cdll.LoadLibrary (os.path.dirname (os.path.realpath (__file__)) + '/lib/libtracer.so')

# Define argument and return types for computeCoefficients
libcoefficients.computeCoefficients.argtypes = [cp.c_double,
                                                cp.c_int, cp.c_int,
                                                np.ctypeslib.ndpointer (cp.c_double,
                                                                        flags = ['C', 'A', 'W', 'O']),
                                                np.ctypeslib.ndpointer (cp.c_double,
                                                                        flags = ['C', 'A', 'W', 'O']),
                                                np.ctypeslib.ndpointer (cp.c_double,
                                                                        flags = ['C', 'A', 'W', 'O']),
                                                np.ctypeslib.ndpointer (cp.c_double,
                                                                        flags = ['C', 'A', 'W', 'O']),
                                                cp.c_int,
                                                np.ctypeslib.ndpointer (cp.c_int,
                                                                        flags = ['C', 'A', 'W', 'O']),
                                                np.ctypeslib.ndpointer (np.complex128, ndim = 2,
                                                                        flags = ['C', 'A', 'W', 'O']),
                                                np.ctypeslib.ndpointer (np.complex128, ndim = 2,
                                                                        flags = ['C', 'A', 'W', 'O']),
                                                np.ctypeslib.ndpointer (np.complex128, ndim = 3,
                                                                        flags = ['C', 'A', 'W', 'O']),
                                                np.ctypeslib.ndpointer (np.complex128, ndim = 3,
                                                                        flags = ['C', 'A', 'W', 'O'])]
libcoefficients.computeCoefficients.restype = None

# Define argument and return types for readModelHeader
libio.readModelHeader.argtypes = [cp.POINTER (cp.c_int),
                                  cp.c_char_p,
                                  cp.c_int,
                                  np.ctypeslib.ndpointer (cp.c_int,
                                                          flags = ['C', 'A', 'W', 'O'])]
libio.readModelHeader.restype = cp.c_int

# Define argument and return types for readModel
libio.readModel.argtypes = [cp.POINTER (cp.c_int),
                            cp.c_int,
                            cp.c_char_p,
                            np.ctypeslib.ndpointer (cp.c_double,
                                                    flags = ['C', 'A', 'W', 'O']),
                            np.ctypeslib.ndpointer (cp.c_double,
                                                    flags = ['C', 'A', 'W', 'O']),
                            np.ctypeslib.ndpointer (cp.c_double,
                                                    flags = ['C', 'A', 'W', 'O']),
                            np.ctypeslib.ndpointer (cp.c_double,
                                                    flags = ['C', 'A', 'W', 'O'])]
libio.readModel.restype = cp.c_int

# Define argument and return types for ttimesAndPaths
libtracer.ttimesAndPaths.argtypes = [cp.c_double,
                                     cp.c_double,
                                     cp.c_char_p,
                                     cp.c_char_p,
                                     cp.c_char_p,
                                     cp.c_int,
                                     cp.c_int,
                                     np.ctypeslib.ndpointer (cp.c_double,
                                                             flags = ['C', 'A', 'W', 'O']),
                                     np.ctypeslib.ndpointer (cp.c_double,
                                                             flags = ['C', 'A', 'W', 'O']),
                                     np.ctypeslib.ndpointer (cp.c_double,
                                                             flags = ['C', 'A', 'W', 'O']),
                                     cp.c_int,
                                     np.ctypeslib.ndpointer (cp.c_int,
                                                             flags = ['C', 'A', 'W', 'O']),
                                     cp.c_int,
                                     np.ctypeslib.ndpointer (cp.c_char, ndim = 2,
                                                             flags = ['C', 'A', 'W', 'O']),
                                     cp.POINTER (cp.c_int),
                                     np.ctypeslib.ndpointer (cp.c_char, ndim = 2,
                                                             flags = ['C', 'A', 'W', 'O']),
                                     np.ctypeslib.ndpointer (cp.c_int,
                                                             flags = ['C', 'A', 'W', 'O']),
                                     np.ctypeslib.ndpointer (cp.c_int,
                                                             flags = ['C', 'A', 'W', 'O']),
                                     np.ctypeslib.ndpointer (cp.c_double, ndim = 2,
                                                             flags = ['C', 'A', 'W', 'O']),
                                     np.ctypeslib.ndpointer (cp.c_double, ndim = 2,
                                                             flags = ['C', 'A', 'W', 'O'])]
libtracer.ttimesAndPaths.restype = cp.c_int

class Constants (object):
  """Loads constants from libio and exposes Python attributes."""
  class _C (cp.Structure):
    """Wrapper for C constants."""
    _fields_ = [('max_string_len',  cp.c_int),
                ('max_phases',      cp.c_int),
                ('max_path_nodes',  cp.c_int),
                ('n_toag',          cp.c_int),
                ('n_coef',          cp.c_int),
                ('n_coes',          cp.c_int),
                ('n_disc',          cp.c_int),
                ('epsilon',         cp.c_double),
                ('earth_radius',    cp.c_double)]

  def __init__ (self, lib):
    """Initializes and fetches constants from C."""
    # Assign library
    self._lib = lib

    # Set function types
    self._lib.getSeisTracerConstants.argtypes = [cp.POINTER (Constants._C)]
    self._lib.getSeisTracerConstants.restype  = None

    # Create structure
    c = Constants._C ()

    # Fetch constants
    self._lib.getSeisTracerConstants (cp.byref (c))

    # Assign constant values
    self.max_string_len  = int (c.max_string_len)
    self.max_phases      = int (c.max_phases)
    self.max_path_nodes  = int (c.max_path_nodes)
    self.n_toag          = int (c.n_toag)
    self.n_coef          = int (c.n_coef)
    self.n_coes          = int (c.n_coes)
    self.n_disc          = int (c.n_disc)
    self.epsilon         = float (c.epsilon)
    self.earth_radius    = float (c.earth_radius)

# Instantiate constants object
_CONSTANTS = Constants (libio)

# Export module-level constants
MAX_STRING_LEN = _CONSTANTS.max_string_len
MAX_PHASES     = _CONSTANTS.max_phases
MAX_PATH_NODES = _CONSTANTS.max_path_nodes
N_TOAG         = _CONSTANTS.n_toag
N_COEF         = _CONSTANTS.n_coef
N_COES         = _CONSTANTS.n_coes
N_DISC         = _CONSTANTS.n_disc
EPSILON        = _CONSTANTS.epsilon
EARTH_RADIUS   = _CONSTANTS.earth_radius
THRESHOLD      = 1e-6

class Phases (object):
  """Creates Phases object given a phases dictionary."""
  def __init__ (self, source_depth, delta, phases_dict, radii):
    """Constructor"""
    # Assign attributes
    self.source_depth  = source_depth
    self.delta         = delta
    self.phases        = phases_dict
    self.radii         = radii

  def __str__ (self):
    """Lists phases included in object."""
    # Format phase list string
    return f"List of phases: {list (self.phases.keys ())}"

  def getInfo (self, phases = 'all', as_string = False):
    """Returns information of a single phase or list of phases."""
    # Build header string
    phases_info = f"#List of phases  (Source depth: {self.source_depth:5.1f} km / Delta:" \
                  f" {self.delta:5.1f}°)  Ttime  Take-off angle   dT/dD     d2T/dD2\n"

    if phases == 'all':
      # Use all phases if not specified
      phases = list (self.phases.keys ())

    # Sort phases by travel time
    phases = sorted (phases,
                     key = lambda i: self.phases[i][0]['ttime']
                     if self.phases[i] else None)

    for phase in phases:

      for arrival in self.phases[phase]:

        ttime    = arrival['ttime']
        take_off = arrival['take-off']
        dT_dD    = arrival['dT_dD']
        d2T_dD2  = arrival['d2T_dD2']

        # Append arrival info
        phases_info += f"{phase:52s} {ttime:11.2f} {take_off:11.3f}" \
                       f" {dT_dD:12.3f} {d2T_dD2:10.3f}\n"

    if as_string:
      # Return as string if requested
      return phases_info[:-1]

    # Print info
    print (phases_info[:-1])

  def getDict (self, phases = 'all'):
    """Returns information of a single phase or list of phases as a dictionary."""
    # Initialize dictionary
    phases_dict = dict ()

    if phases == 'all':
      # Use all phases if not specified
      phases = list (self.phases.keys ())

    for phase in phases:

      if phase not in phases_dict:
        # Initialize list for phase
        phases_dict[phase] = list ()

      for arrival in self.phases[phase]:

        ttime     = arrival['ttime']
        take_off  = arrival['take-off']
        dT_dD     = arrival['dT_dD']
        d2T_dD2   = arrival['d2T_dD2']
        amplitude = arrival['amplitude']

        # Append arrival dictionary
        phases_dict[phase] += [{'ttime' : ttime,
                                'take-off' : take_off,
                                'dT_dD' : dT_dD,
                                'd2T_dD2' : d2T_dD2,
                                'amplitude' : amplitude}]

    return phases_dict

  def getPath (self, phases = 'all'):
    """Returns ray paths of a single phase or list of phases."""
    # Initialize paths dictionary
    phases_paths = dict ()

    if phases == 'all':
      # Use all phases if not specified
      phases = list (self.phases.keys ())

    for phase in phases:

      if phase not in phases_paths:
        # Initialize list for phase
        phases_paths[phase] = list ()

      for arrival in self.phases[phase]:

        # Append path dictionary
        phases_paths[phase] += [{'radii' : arrival['radii'],
                                 'deltas' : arrival['deltas'],
                                 'branches' : arrival['branches']}]

    return phases_paths

  def plotPath (self, phases = 'all', only_first_arrival = False, save_figure = False):
    """Plot ray paths with branch colors using LineCollection."""
    if phases == 'all':
      # Use all phases if not specified
      phases = list (self.phases.keys ())

    # Define branch colors
    branch_colors = {
      'P': 'red',   'p': 'lightcoral',
      'S': 'blue',  's': 'lightblue',
      'K': 'orange',
      'I': 'purple',
      'J': 'green'
    }

    if save_figure:
      # Use non-interactive backend
      matplotlib.use ('Agg')
    
    else:
      # Use default backend
      matplotlib.use (DEFAULT_BACKEND)

    # Create figure and axis
    fig = plt.figure (figsize = (7, 7), dpi = 140,
                      num = 1, clear = True)
    ax = fig.add_subplot (111)

    # Generate theta for circles
    theta = np.linspace (0, 2 * np.pi, 1000)
    
    cos_theta, sin_theta = np.cos (theta), np.sin (theta)
    
    for radius in self.radii:
      # Plot discontinuity circles
      ax.plot (radius * cos_theta,
               radius * sin_theta,
               color = 'black',
               linewidth = 0.2)

    # Initialize segments and colors
    all_segments = list ()
    all_colors = list ()

    for phase in phases:

      for arrival in self.phases[phase]:
      
        Rs = arrival['radii']
        Ds = arrival['deltas']
        Br = arrival['branches']

        # Compute Cartesian coordinates
        cos_Ds, sin_Ds = np.cos(Ds), np.sin(Ds)
        Xs, Ys = Rs * cos_Ds, Rs * sin_Ds

        # Create segment points
        p0 = np.column_stack ((Xs[:-1], Ys[:-1]))
        p1 = np.column_stack ((Xs[1:],  Ys[1:]))
        
        segs = np.stack ((p0, p1), axis = 1)
        all_segments.append (segs)

        # Assign segment colors
        colors = [branch_colors.get (branch, 'gray') for branch in Br[:-1]]
        all_colors.extend (colors)

        # Get source and receiver positions
        xs, ys = Xs[0],  Ys[0]
        xe, ye = Xs[-1], Ys[-1]

        if only_first_arrival:
          # Stop after first arrival
          break

    if all_segments:
      
      # Concatenate segments
      all_segments = np.concatenate (all_segments, axis = 0)
      # Create line collection
      lc = LineCollection (all_segments,
                           linewidths = 0.8,
                           colors = all_colors)
      ax.add_collection (lc)

    # Plot source
    ax.scatter (xs, ys, marker = '*',
                color = 'yellow', edgecolor = 'black',
                linewidth = 0.5, s = 150, zorder = 100)
    # Plot receiver
    ax.scatter (xe, ye, marker = (3, 0, 0),
                color = 'red', edgecolor = 'black',
                linewidth = 0.5, s = 150, zorder = 100)

    # Create legend handles
    legend_handles = [
      Line2D ([0], [0], color = color, lw = 2, label = branch)
      for branch, color in branch_colors.items ()
    ]
    # Add legend
    ax.legend (handles = legend_handles,
               loc = 'upper right',
               frameon = True,
               fontsize = 8)

    # Set plot limits
    plt.xlim (-6500, 6500)
    plt.ylim (-6500, 6500)
    plt.xticks ([])
    plt.yticks ([])

    # Adjust subplot
    plt.subplots_adjust (left = 0.0, bottom = 0.0,
                         top = 1.0, right = 1.0)
    plt.axis ('off')

    if save_figure:
      # Save figure
      plt.savefig (f'ray_paths_{self.source_depth:.0f}_{self.delta:.0f}.png',
                   bbox_inches = 'tight')

    else:
      # Show plot
      plt.show ()

    # Close plot
    plt.close ()

class Tracer (object):
  """Creates Tracer object given a velocity model. Default is 'IASP91'."""
  def __init__ (self, model = 'iasp91'):
    """Constructor"""
    # Assign model name
    self.model = model

    # Find available models
    paths = glob (os.path.dirname (os.path.realpath (__file__)) + '/models/*.model')
    models = [Path(p).stem for p in paths]

    if self.model not in models:
      
      sys.exit (f"Error: model '{self.model}' is not yet implemented! Available models are: {models}")

    # Set model and tables paths
    self._model_path  = os.path.dirname (os.path.realpath (__file__)) + f'/models/{self.model}.model'
    self._tables_path = os.path.dirname (os.path.realpath (__file__)) + f'/tables/{self.model}'

    # Initialize variables
    self._nl      = cp.c_int ()
    self._N_TOAG  = N_TOAG
    self._N_COEF  = N_COEF
    self._N_COES  = N_COES
    self._N_DISC  = N_DISC
    self._nd      = self._N_DISC + 2
    self._indices = np.empty (self._nd, dtype = np.int32)

    # Read model header
    libio.readModelHeader (cp.byref (self._nl),
                           cp.c_char_p (self._model_path.encode ('utf-8')),
                           cp.c_int (self._nd),
                           self._indices)

    # Initialize model arrays
    self._ni   = cp.c_int ()
    self._rn   = np.empty (self._nl.value, dtype = np.float64)
    self._rhon = np.empty (self._nl.value, dtype = np.float64)
    self._vpn  = np.empty (self._nl.value, dtype = np.float64)
    self._vsn  = np.empty (self._nl.value, dtype = np.float64)

    # Read model data
    return_value = libio.readModel (cp.byref (self._ni),
                                    self._nl.value,
                                    cp.c_char_p (self._model_path.encode ('utf-8')),
                                    self._rn,
                                    self._rhon,
                                    self._vpn,
                                    self._vsn)

    if return_value == 1:

      sys.exit (f"Error: could not open file '{self._model_path}'!")

    # Define coefficient keys
    self._keysF = {'SHfSH' : 0, 'PfP' : 1, 'PfSV' : 2, 'SVfP' : 3, 'SVfSV' : 4}
    self._keysS = {'SHdSHd' : 0, 'SHdSHu' : 1, 'SHuSHd' : 2, 'SHuSHu' : 3,
                   'PdPd' : 4, 'PdPu' : 5, 'PuPd' : 6, 'PuPu' : 7,
                   'PdSVd' : 8, 'PdSVu' : 9, 'PuSVd' : 10, 'PuSVu' : 11,
                   'SVdSVd' : 12, 'SVdSVu' : 13, 'SVuSVd' : 14, 'SVuSVu' : 15,
                   'SVdPd' : 16, 'SVdPu' : 17, 'SVuPd' : 18, 'SVuPu' : 19}

    # Define discontinuity map
    self._map = {'surface' : 0, 'moho' : 1,
                 'd410' : 2, 'd660' : 3,
                 'cmb' : 4, 'icb' : 5}

  def __str__ (self):
    """Prints velocity model name and parameters."""
    # Format model string
    return f"Model object: {self.model}\n\nR (km): {self._rn}\n\n" \
           f"Rho (g/cm3): {self._rhon}\n\nVp (km/s): {self._vpn}\n\nVs (km/s): {self._vsn}"

  def _findDeltasForPhases (self, source_depth = 0, requested_phases = ['P', 'S']):
    """Collects all unique deltas from the table for each requested phase."""
    # Copy original phases
    original_phases  = copy.deepcopy (requested_phases)
    difracted_phases = list ()
  
    if 'Pdiff' in requested_phases:
      # Add Pdiff if requested
      difracted_phases.append ('Pdiff')
      requested_phases.remove ('Pdiff')

      if 'P' not in requested_phases:

        requested_phases.append ('P')

    if 'Sdiff' in requested_phases:
      # Add Sdiff if requested
      difracted_phases.append ('Sdiff')
      requested_phases.remove ('Sdiff')

      if 'S' not in requested_phases:

        requested_phases.append ('S')

    # Round source depth
    dp = int (round (source_depth))

    # Set paths
    phases_path = f"{self._tables_path}/phases_{dp}.list"
    table_path = f"{self._tables_path}/table_{dp}.bin"

    # Open phases file
    with open (phases_path, 'r') as f:

      lines = f.readlines ()

    # Parse number of phases
    nP = int (lines[0].strip ())
    phase_data = lines[1:]
    
    if len (phase_data) != nP:
      
      sys.exit ('Mismatch in number of phases.')

    # Initialize lists
    Nl_list = list ()
    phase_list = list ()
    
    for line in phase_data:
    
      parts = line.strip ().split (maxsplit = 1)
    
      if len (parts) != 2:
    
        raise sys.exit ('Invalid phase line.')
      
      # Parse phase data
      Nl, phase = int (parts[0]), parts[1]
      Nl_list.append (Nl)
      phase_list.append (phase)

    # Find requested indices
    requested_indices = []
    
    for i, phase in enumerate (phase_list):
      
      if phase in requested_phases:
        
        requested_indices.append (i)

    if not requested_indices:
      
      return {}

    # Open table file
    with open (table_path, 'rb') as f:
      
      data = f.read ()

    offset = 0
    # Unpack number of blocks
    N = struct.unpack_from ('i', data, offset)[0]
    offset += 4

    if len (data) - 4 != N * 16:
      
      sys.exit ('Incomplete or invalid table data.')

    # Parse blocks
    blocks = list ()
    
    for i in range (N):
  
      block = struct.unpack_from ('ffff', data, offset)
      np_ = int (block[0])
      index = int (block[1])
      delta = float (block[2])
      psph = float (block[3])
      blocks.append ((np_, index, delta, psph))
      
      offset += 16

    # Compute block starts
    block_starts = [0]
    
    for nl in Nl_list[:-1]:
    
      block_starts.append (block_starts[-1] + nl)

    # Initialize phase deltas
    phase_deltas = {phase: set () for phase in requested_phases}

    for idx in requested_indices:
 
      phase = phase_list[idx]
      start = block_starts[idx]
      end = start + Nl_list[idx]

      for b in range (start, end):

        delta = blocks[b][2]
        psph = blocks[b][3]
        # Add delta and psph
        phase_deltas[phase].add ((delta, psph))

    # Sort phase deltas
    phase_deltas = {phase: np.array (sorted (list (deltas),
                                     key = lambda x: (x[1], x[0])))
                                     for phase, deltas in phase_deltas.items ()}

    if 'Pdiff' in difracted_phases:

      # Set Pdiff deltas
      delta_min = phase_deltas['P'][0, 0]
      delta_max = 150.0
      dT_dD     = phase_deltas['P'][0, 1]

      phase_deltas['Pdiff']       = np.empty ((2, 2), dtype = np.float64)
      phase_deltas['Pdiff'][0, 0] = delta_max
      phase_deltas['Pdiff'][1, 0] = delta_min
      phase_deltas['Pdiff'][0, 1] = dT_dD
      phase_deltas['Pdiff'][1, 1] = dT_dD

      if 'P' not in original_phases:
        # Remove P if not original
        del phase_deltas['P']

    if 'Sdiff' in difracted_phases:

      # Set Sdiff deltas
      delta_min = phase_deltas['S'][0, 0]
      delta_max = 150.0
      dT_dD     = phase_deltas['S'][0, 1]

      phase_deltas['Sdiff']       = np.empty ((2, 2), dtype = np.float64)
      phase_deltas['Sdiff'][0, 0] = delta_max
      phase_deltas['Sdiff'][1, 0] = delta_min
      phase_deltas['Sdiff'][0, 1] = dT_dD
      phase_deltas['Sdiff'][1, 1] = dT_dD

      if 'S' not in original_phases:
        # Remove S if not original
        del phase_deltas['S']

    return phase_deltas

  def _mergePhases (self,
                    n_phases_1, phases_list_1, Ni_1, Rs_1, Ds_1,
                    n_phases_2, phases_list_2, Ni_2, Rs_2, Ds_2):
    """Merges two ordered sets of seismic phases."""      
    # Calculate total phases
    n_phases = n_phases_1 + n_phases_2

    # Initialize merged arrays
    phases_list = np.zeros ((n_phases, MAX_STRING_LEN), dtype = 'S1')

    Ni = np.zeros (n_phases, dtype = np.int32)

    Rs = np.zeros ((n_phases, MAX_PATH_NODES))
    Ds = np.zeros ((n_phases, MAX_PATH_NODES))

    # Extract times
    times_1 = [phases_list_1[i].tobytes ().split ()[1] for i in range (n_phases_1)]
    times_2 = [phases_list_2[i].tobytes ().split ()[1] for i in range (n_phases_2)]

    # Merge loop variables
    i = 0; j = 0; k = 0

    while i < n_phases_1 and j < n_phases_2:

      if times_1[i] <= times_2[j]:

        # Copy from first set
        phases_list[k] = phases_list_1[i]
            
        Ni[k] = Ni_1[i]
            
        Rs[k] = Rs_1[i]
        Ds[k] = Ds_1[i]
            
        i += 1

      else:
            
        # Copy from second set
        phases_list[k] = phases_list_2[j]
            
        Ni[k] = Ni_2[j]
            
        Rs[k] = Rs_2[j]
        Ds[k] = Ds_2[j]
            
        j += 1

      k += 1

    while i < n_phases_1:
          
      # Copy remaining from first
      phases_list[k] = phases_list_1[i]
          
      Ni[k] = Ni_1[i]

      Rs[k] = Rs_1[i]
      Ds[k] = Ds_1[i]
          
      i += 1; k += 1

    while j < n_phases_2:
          
      # Copy remaining from second
      phases_list[k] = phases_list_2[j]
          
      Ni[k] = Ni_2[j]

      Rs[k] = Rs_2[j]
      Ds[k] = Ds_2[j]
          
      j += 1; k += 1

    return n_phases, phases_list, Ni, Rs, Ds

  def _mergeMinorAndMajor (self,
                           n_phases_m, phases_list_m, Ni_m, Rs_m, Ds_m,
                           n_phases_M, phases_list_M, Ni_M, Rs_M, Ds_M):
    """Merges minor-arc and major-arc seismic phases."""      
    # Extract names
    names_m = [phases_list_m[i].tobytes ().split ()[0] for i in range (n_phases_m)]
    names_M = [phases_list_M[i].tobytes ().split ()[0] for i in range (n_phases_M)]

    # Find duplicates
    duplicates = set (names_m).intersection (names_M)

    w = 48

    for i in range (n_phases_m):
          
      name = names_m[i]
          
      if name in duplicates:
          
        # Append suffix for minor
        name += b'_m'
        phases_list_m[i, :w] = np.frombuffer (name[:w].ljust (w, b' '),
                                              dtype = 'S1')

    for i in range (n_phases_M):
          
      name = names_M[i]
          
      if name in duplicates:
          
        # Append suffix for major
        name += b'_M'
        phases_list_M[i, :w] = np.frombuffer (name[:w].ljust (w, b' '),
                                              dtype = 'S1')

    # Merge phases
    return self._mergePhases (n_phases_m, phases_list_m, Ni_m, Rs_m, Ds_m,
                              n_phases_M, phases_list_M, Ni_M, Rs_M, Ds_M)
    
  def _boundariesFromModel (self):
    """Major discontinuities (depths, km) derived from model radii."""
    # Compute radii and depths
    radii  = self._rn[self._indices].astype (float)
    R0     = float (radii[0])
    depths = R0 - radii

    # Return boundaries dictionary
    return {'surface' : float (depths[0]),
            'moho'    : float (depths[1]),
            'd410'    : float (depths[2]),
            'd660'    : float (depths[3]),
            'cmb'     : float (depths[4]),
            'icb'     : float (depths[5]),
            'center'  : float (depths[6])}

  def _edges (self, boundaries):
    """Ordered list of shell edges (top -> bottom)."""
    # Return list of edges
    return [boundaries['surface'],
            boundaries['moho'],
            boundaries['d410'],
            boundaries['d660'],
            boundaries['cmb'],
            boundaries['icb'],
            boundaries['center']]

  def _isBoundaryPoint (self, Rs):
    """True if the ray path is on a shell edge."""
    # Get boundary radii
    boundary_radii = self._rn[self._indices[:-1]].astype (float)

    # Check closeness
    is_boundary_point = np.isclose (Rs[:, None],
                                    boundary_radii,
                                    atol = EPSILON,
                                    rtol = 0.0).any (axis = 1)

    if len (Rs) >= 1:
      # Exclude endpoints
      is_boundary_point[0]  = False
      is_boundary_point[-1] = False

    return is_boundary_point

  def _shellBoundsFromDepthDir (self, depth, going_up, boundaries):
    """(top, bottom) of the shell that contains depth.
       If exactly on an interface, tie-break with direction."""
    # Get edges
    edges = self._edges (boundaries)

    for k, edge in enumerate (edges):
      
      if abs (depth - edge) <= EPSILON:
      
        if going_up:
      
          # Adjust for upward
          k1 = max (0, k - 1)
          
          return (edges[k1], edges[k])
      
        else:
      
          # Adjust for downward
          k2 = min (len (edges) - 1, k + 1)
          
          return (edges[k], edges[k2])

    for top, bottom in zip (edges[:-1], edges[1:]):
      
      if top <= depth <= bottom:
      
        return top, bottom

    # Default bounds
    return (edges[0], edges[1]) if going_up else (edges[-2], edges[-1])

  def _nextBoundary (self, depth, going_up, boundaries):
    """Nearest boundary above/below the given depth."""
    # Get edges
    edges = self._edges (boundaries)

    if going_up:
      
      # Find shallower boundaries
      shallower = [edge for edge in edges if edge < depth]
    
      return max (shallower) if shallower else boundaries['surface']

    # Find deeper boundaries
    deeper = [edge for edge in edges if edge > depth]
    
    return min (deeper) if deeper else boundaries['center']

  def _decodeCodeToLegs (self, code, source_depth):
    """Turn code into depth-monotone legs using Pu/Su/Pd/Sd mapping."""
    # Get boundaries
    boundaries = self._boundariesFromModel ()

    # Initialize depth and legs
    depth = float (source_depth)
    legs  = list ()

    lower_case = None

    for i, c in enumerate (code):

      # Determine branch and direction
      branch = 'P' if c in ['A', 'C'] else 'S'
      going_up = True if c in ['A', 'B'] else False

      if i == 0:
        
        lower_case = going_up

      if lower_case and branch in ['P', 'S']:
        
        # Lowercase for upward
        branch = branch.lower ()

      # Get shell bounds
      top, bottom = self._shellBoundsFromDepthDir (depth,
                                                   going_up,
                                                   boundaries)
      # Get next boundary
      boundary_depth = self._nextBoundary (depth,
                                           going_up,
                                           boundaries)

      # Append leg dict
      legs.append (dict (branch = branch,
                         up     = going_up,
                         top    = top,
                         bottom = bottom))

      if lower_case and going_up and i + 1 < len (code) and \
        not (True if code[i + 1] in ['A', 'B'] else False):
        
        # Reset lower case
        lower_case = False

      # Update depth
      depth = boundary_depth

    return legs, boundaries

  def _classifyBranchByDepth (self, branch, depth, boundaries):
    """Return the phase label for a given branch and depth.
       Labels are:
         P/p : P-wave in crust or mantle (lowercase only for first upward leg)
         S/s : S-wave in crust or mantle (lowercase only for first upward leg)
         K   : P-wave in outer core
         I   : P-wave in inner core
         J   : S-wave in inner core"""
    # Get boundaries
    cmb, icb = boundaries['cmb'], boundaries['icb']

    if branch.upper () == 'P':
      
      # Classify P branch
      return 'K' if cmb <= depth < icb else \
             'I' if depth >= icb else branch

    if branch.upper () == 'S':
      
      # Classify S branch
      return 'J' if depth >= icb else branch

    return branch

  def _branchesAlongPath (self, code, source_depth, Rs, phase):
    """Assign a wave-type label (P, p, S, s, K, I, J) to each path node."""
    # Decode legs
    legs, boundaries = self._decodeCodeToLegs (code, source_depth)

    # Compute depths
    depths = boundaries['center'] - Rs
    is_boundary_point = self._isBoundaryPoint (Rs)

    # Initialize leg indices
    N = len (Rs); i = 0; li = 0; L = np.zeros (N, dtype = np.int32)

    while i < N:

      L[i] = li

      if i > 0 and is_boundary_point[i - 1]:
        
        L[i - 1] = li

      if is_boundary_point[i]:
        
        # Increment for boundary
        li += 1; i += 1
      
      elif i > 0 and depths[i - 1] < depths[i] \
                 and depths[i + 1] < depths[i]:
        
        # Increment for turning point
        li += 1

      i += 1

    # Initialize branches array
    branches = np.empty (N, dtype = '<U1')

    for i in range (N):

      # Classify each node
      branches[i] = self._classifyBranchByDepth (legs[L[i]]['branch'],
                                                 depths[i],
                                                 boundaries)

    return branches

  def _classifyDiscontinuities (self, radii, branches):
    """Classify all discontinuities crossed by the ray and return the interaction type.
       Surface reflections:   SHfSH, PfP, PfSV, SVfP, SVfSV
       Internal reflections/transmissions:
         SHdSHd, SHdSHu, SHuSHd, SHuSHu,
         PdPd, PdPu, PuPd, PuPu,
         PdSVd, PdSVu, PuSVd, PuSVu,
         SVdSVd, SVdSVu, SVuSVd, SVuSVu,
         SVdPd, SVdPu, SVuPd, SVuPu
       
       A branch is treated as SH when the entire ray path consists only of
       shear-wave legs (s, S, J). Otherwise SV is assumed."""
    # Get boundaries
    boundaries = self._boundariesFromModel ()

    # Check if SH wave
    is_sh = np.all (np.isin (branches, ['s', 'S', 'J']))
    # Find reflection indices
    indices = np.flatnonzero (np.diff (radii) == 0)[1:]

    # Initialize interactions
    interactions = list ()

    for index in indices:

      # Compute depth and find disc
      depth = EARTH_RADIUS - radii[index]
      disc  = min (boundaries, key = lambda k: abs (boundaries[k] - depth))

      key = ''

      if np.isin (branches[index - 1], ['s', 'S', 'J']):

        # Determine incoming type
        key += 'SH' if is_sh else 'SV'

      else:

        key += 'P'

      if abs (depth) < THRESHOLD:

        # Surface reflection
        key += 'f'

      elif radii[index - 1] < radii[index]:

        key += 'u'

      else:

        key += 'd'

      if np.isin (branches[index + 1], ['s', 'S', 'J']):

        # Determine outgoing type
        key += 'SH' if is_sh else 'SV'

      else:

        key += 'P'

      if abs (depth) < THRESHOLD:

        pass

      elif radii[index + 1] < radii[index + 2]:

        key += 'u'

      else:

        key += 'd'

      # Append interaction
      interactions.append ((disc, key))
    
    return interactions

  def plotModel (self):
    """Plots the loaded Earth model."""
    # Compute depth
    depth = EARTH_RADIUS - self._rn
    disc_depths = EARTH_RADIUS - self._rn[self._indices[1:-1]]

    # Define labels
    labels = ['Moho', 'd410', 'd660', 'CMB', 'ICB']

    # Create figure
    plt.figure (figsize = (5, 6), dpi = 140)

    # Plot velocities and density
    plt.plot (self._vpn, depth, label = r'$\alpha$ [km/s]', color = 'red')
    plt.plot (self._vsn, depth, label = r'$\beta$ [km/s]', color = 'blue')
    plt.plot (self._rhon, depth, label = r'$\rho$ [g/cm³]', color = 'green')

    for d_depth, label in zip (disc_depths, labels):
      
      # Plot discontinuities
      plt.axhline (d_depth, color = 'black', linestyle = 'dashed',
                   linewidth = 1)
      plt.text (plt.gca ().get_xlim ()[1] * 1.02, d_depth, label,
                fontsize = 12, color = 'black')

    # Set labels and limits
    plt.xlabel ('Velocity ([km/s] / Density [g/cm³]', fontsize = 12)
    plt.ylabel ('Depth [km]', fontsize = 12)
    plt.ylim (EARTH_RADIUS, 0)
    plt.tick_params (axis = 'x', labelsize = 12)
    plt.tick_params (axis = 'y', labelsize = 12)
    plt.grid (linestyle = 'dashed', linewidth = '0.5')
    plt.legend ()
    plt.title (f'{self.model.upper ()} Earth Model', fontsize = 14)
    plt.tight_layout ()
    plt.show ()

  def computeAmplitudeCoefficients (self):
    """Computes amplitude coefficients at all the internal discontinuities."""
    # Get boundaries
    boundaries = self._boundariesFromModel ()

    # Initialize coefficients dict
    coefficients = dict ()

    for disc in boundaries.keys ():

      if disc == 'center':
        # Skip center
        continue

      # Set source depth
      source_depth = boundaries[disc]

      # Initialize coefficient arrays
      cf_FP = np.zeros ((self._N_TOAG, self._N_COEF), dtype = np.complex128)
      cf_FS = np.zeros ((self._N_TOAG, self._N_COEF), dtype = np.complex128)
      cf_SP = np.zeros ((self._N_TOAG, self._N_DISC, self._N_COES), dtype = np.complex128)
      cf_SS = np.zeros ((self._N_TOAG, self._N_DISC, self._N_COES), dtype = np.complex128)

      coefficients[disc] = dict ()

      if disc == 'surface':

        # Adjust for surface
        source_depth += THRESHOLD

        # Compute coefficients
        libcoefficients.computeCoefficients (source_depth, self._ni.value, self._nl.value,
                                             self._rn, self._rhon, self._vpn, self._vsn,
                                             self._nd, self._indices,
                                             cf_FP, cf_FS, cf_SP, cf_SS)

        for key in self._keysF.keys ():

          key_index = self._keysF[key]

          if key[0] == 'P':

            # Copy P coefficients
            coefficients[disc][key] = cf_FP[:, key_index].copy ()
          
          else:
          
            # Copy S coefficients
            coefficients[disc][key] = cf_FS[:, key_index].copy ()

      else:

        # Get disc index
        disc_index = self._map[disc]

        # Separate keys
        keys1 = [key for key in self._keysS.keys () if key[1] == 'd' or  key[2] == 'd']
        keys2 = [key for key in self._keysS.keys () if key[1] == 'u' or  key[2] == 'u']

        # Adjust depths
        source_depth1 = source_depth - THRESHOLD
        source_depth2 = source_depth + THRESHOLD

        # Compute for first depth
        libcoefficients.computeCoefficients (source_depth1, self._ni.value, self._nl.value,
                                             self._rn, self._rhon, self._vpn, self._vsn,
                                             self._nd, self._indices,
                                             cf_FP, cf_FS, cf_SP, cf_SS)

        for key in keys1:

          key_index = self._keysS[key]

          if key[0] == 'P':

            # Copy P coefficients
            cf = cf_SP[:, disc_index - 1, key_index].copy ()
          
          else:
          
            # Copy S coefficients
            cf = cf_SS[:, disc_index - 1, key_index].copy ()

          coefficients[disc][key] = cf 

        # Compute for second depth
        libcoefficients.computeCoefficients (source_depth2, self._ni.value, self._nl.value,
                                             self._rn, self._rhon, self._vpn, self._vsn,
                                             self._nd, self._indices,
                                             cf_FP, cf_FS, cf_SP, cf_SS)

        for key in keys2:

          key_index = self._keysS[key]

          if key[0] == 'P':

            # Copy P coefficients
            cf = cf_SP[:, disc_index - 1, key_index].copy ()
          
          else:
          
            # Copy S coefficients
            cf = cf_SS[:, disc_index - 1, key_index].copy ()

          coefficients[disc][key] = cf 

    return coefficients

  def plotAmplitudeCoefficients (self, discontinuities = 'all'):
    """Plots amplitude coefficients at all the internal discontinuities."""
    # Get boundaries
    boundaries = list (self._boundariesFromModel ().keys ())[:-1]

    if discontinuities == 'all':
      
      # Use all discontinuities
      disc_keys = boundaries[::]

    else:

      if type (discontinuities) is str:
        # Convert to list
        discontinuities = [discontinuities]

      # Validate discontinuities
      disc_keys = list ()

      for disc in discontinuities:

        if disc in boundaries:

          disc_keys += [disc]

        else:

          sys.exit (f"Error: discontinuities must be on the list: {boundaries}.")

    # Compute coefficients
    coefficients = self.computeAmplitudeCoefficients ()

    # Generate theta and sin_theta
    theta = np.linspace (0, 90, N_TOAG)
    sin_theta = np.sin (np.radians (theta))

    # Define key groups
    keys_Pf  = ['PfP', 'PfSV']
    keys_SVf = ['SVfSV', 'SVfP']
    keys_SHf = ['SHfSH']

    keys_Pd  = ['PdPd', 'PdSVd', 'PdPu', 'PdSVu']
    keys_Pu  = ['PuPu', 'PuSVu', 'PuPd', 'PuSVd']
    keys_SVd = ['SVdSVd', 'SVdPd', 'SVdSVu', 'SVdPu']
    keys_SVu = ['SVuSVu', 'SVuPu', 'SVuSVd', 'SVuPd']
    keys_SHd = ['SHdSHd', 'SHdSHu']
    keys_SHu = ['SHuSHu', 'SHuSHd']

    for disc in disc_keys:

      # Get index
      index = self._indices[self._map[disc]]

      if disc == 'surface':
      
        # Surface keys
        all_keys = [keys_Pf, keys_SVf, keys_SHf]
      
      else:

        # Internal keys
        all_keys = [keys_Pd, keys_Pu,
                    keys_SVd, keys_SVu,
                    keys_SHd, keys_SHu]

      for keys in all_keys:

        # Create subplots
        fig, axs = plt.subplots (1, 2, dpi = 140)

        # Initialize flux
        total_energy_flux = np.zeros_like (sin_theta)

        for key in keys:

          # Get coefficient
          cf = coefficients[disc][key]

          if disc == 'surface':

            # Get velocities
            alpha = self._vpn[index]
            beta  = self._vsn[index] 
          
            # Compute p
            p = sin_theta / alpha if key[0] == 'P' else sin_theta / beta

            # Compute ci, cj
            ci = 1.0 - (alpha * p) ** 2
            cj = 1.0 - (beta  * p) ** 2

            # Compute cosines
            cosi = np.real (np.sqrt (ci.astype (complex)))
            cosj = np.real (np.sqrt (cj.astype (complex)))

            # Compute k values
            ki = alpha * cosi
            kj = beta  * cosj

            # Compute amplitude and energy
            cf_A = np.abs (cf)
            cf_E = cf_A ** 2

            # Compute ratios
            ra = np.divide (kj, ki, where = (ki > EPSILON),
                                    out = np.zeros_like (ki))
            rb = np.divide (ki, kj, where = (kj > EPSILON),
                                    out = np.zeros_like (kj))

            if   key == 'PfSV': cf_E *= ra
            elif key == 'SVfP': cf_E *= rb

          else:

            # Get layer velocities and densities
            alpha1 = self._vpn[index - 1]
            alpha2 = self._vpn[index]
            beta1  = self._vsn[index - 1]
            beta2  = self._vsn[index]
            rho1   = self._rhon[index - 1]
            rho2   = self._rhon[index]

            # Set kc
            kc = 1.0

            if ((key[:3] == 'SVd' or key[:3] == 'SHd')
                and beta1 < EPSILON) or \
               ((key[:3] == 'SVu' or key[:3] == 'SHu')
                and beta2 < EPSILON) or \
               ((key[-3:] == 'SVd' or key[-3:] == 'SHd')
                and beta2 < EPSILON) or \
               ((key[-3:] == 'SVu' or key[-3:] == 'SHu')
                and beta1 < EPSILON):

              kc = 0.0

            if key[1] == 'd' or key[2] == 'd':
              
              # Set alpha, beta for down
              alpha = alpha1
              beta  = beta1

            else:
            
              # Set alpha, beta for up
              alpha = alpha2
              beta  = beta2

            if beta < EPSILON:
              
              # Avoid zero beta
              beta = EPSILON

            # Compute p
            p = sin_theta / alpha if key[0] == 'P' else sin_theta / beta

            # Compute ci, cj for layers
            ci1 = 1.0 - (alpha1 * p) ** 2
            ci2 = 1.0 - (alpha2 * p) ** 2
            cj1 = 1.0 - (beta1  * p) ** 2
            cj2 = 1.0 - (beta2  * p) ** 2

            # Compute cosines
            cosi1 = np.real (np.sqrt (ci1.astype (complex)))
            cosi2 = np.real (np.sqrt (ci2.astype (complex)))
            cosj1 = np.real (np.sqrt (cj1.astype (complex)))
            cosj2 = np.real (np.sqrt (cj2.astype (complex)))

            # Compute k and kr values
            ki1 = alpha1 * cosi1; kri1 = rho1 * ki1
            ki2 = alpha2 * cosi2; kri2 = rho2 * ki2
            kj1 = beta1  * cosj1; krj1 = rho1 * kj1
            kj2 = beta2  * cosj2; krj2 = rho2 * kj2

            # Compute amplitude and energy
            cf_A = kc * np.abs (cf)
            cf_E = cf_A ** 2

            # Compute ratios
            ra = np.divide (kri2, kri1, where = (kri1 > EPSILON),
                                        out = np.zeros_like (kri2))
            rb = np.divide (kj1,  ki1,  where = (ki1  > EPSILON),
                                        out = np.zeros_like (kj1))
            rc = np.divide (krj2, kri1, where = (kri1 > EPSILON),
                                        out = np.zeros_like (krj2))
            rd = np.divide (kri1, kri2, where = (kri2 > EPSILON),
                                        out = np.zeros_like (kri1))
            re = np.divide (kj2,  ki2,  where = (ki2  > EPSILON),
                                        out = np.zeros_like (kj2))
            rf = np.divide (krj1, kri2, where = (kri2 > EPSILON),
                                        out = np.zeros_like (krj1))
            rg = np.divide (ki1,  kj1,  where = (kj1  > EPSILON),
                                        out = np.zeros_like (ki1))
            rh = np.divide (kri2, krj1, where = (krj1 > EPSILON),
                                        out = np.zeros_like (kri2))
            ri = np.divide (krj2, krj1, where = (krj1 > EPSILON),
                                        out = np.zeros_like (krj2))
            rj = np.divide (kri1, krj2, where = (krj2 > EPSILON),
                                        out = np.zeros_like (kri1))
            rk = np.divide (ki2,  kj2,  where = (kj2  > EPSILON),
                                        out = np.zeros_like (ki2))
            rl = np.divide (krj1, krj2, where = (krj2 > EPSILON),
                                        out = np.zeros_like (krj1))

            # Adjust energy flux
            if   key == 'PdPd':   cf_E *= ra
            elif key == 'PdSVu':  cf_E *= rb
            elif key == 'PdSVd':  cf_E *= rc
            elif key == 'PuPu':   cf_E *= rd
            elif key == 'PuSVd':  cf_E *= re
            elif key == 'PuSVu':  cf_E *= rf
            elif key == 'SVdPu':  cf_E *= rg
            elif key == 'SVdPd':  cf_E *= rh
            elif key == 'SVdSVd' or \
                 key == 'SHdSHd': cf_E *= ri
            elif key == 'SVuPu':  cf_E *= rj
            elif key == 'SVuPd':  cf_E *= rk
            elif key == 'SVuSVu' or \
                 key == 'SHuSHu': cf_E *= rl

          # Add to total flux
          total_energy_flux += cf_E

          # Plot amplitude
          axs[0].plot (theta, cf_A)
          # Plot energy
          axs[1].plot (theta, cf_E, label = f'{key}')

        # Plot total flux
        axs[1].plot (theta, total_energy_flux, color = 'black',
                     label = 'Total')
        # Set titles
        axs[0].set_title ('Amplitude Coefficients (Norm)',
                          fontsize = 12, y = 1.02)
        axs[1].set_title ('Energy Flux',
                          fontsize = 12, y = 1.02)
        # Set ticks
        axs[0].tick_params (axis = 'x', labelsize = 12)
        axs[0].tick_params (axis = 'y', labelsize = 12)
        axs[1].tick_params (axis = 'x', labelsize = 12)
        axs[1].tick_params (axis = 'y', labelsize = 12)
        # Set limits
        axs[0].set_xlim (0, 90)
        axs[1].set_xlim (0, 90)
        # Set labels
        axs[0].set_xlabel (r'Angle of incidence [$\circ$]',
                           fontsize = 12)
        axs[1].set_xlabel (r'Angle of incidence [$\circ$]',
                           fontsize = 12)
        # Add legend
        axs[1].legend (fontsize = 12)
        # Add grids
        axs[0].grid (linestyle = 'dashed', linewidth = '0.5')
        axs[1].grid (linestyle = 'dashed', linewidth = '0.5')
        # Set suptitle
        plt.suptitle (disc.upper (), fontsize = 14, y = 0.96)
        # Adjust layout
        fig.tight_layout (pad = 1.4)

      # Show plots
      plt.show ()

  def _ttimesAndPaths (self, source_depth = 0, delta = None, phases = 'all', arc = 'minor',
                       min_loop = 0, max_loop = 0, compute_amplitudes = True):
    """Computes seismic waves travel times and ray paths given a source depth and
       source-receiver distance."""
    # Validate delta
    if delta == None or delta < 0.0 or delta > 180.0:

      sys.exit ("Error: 'delta' (source-receiver distance) should be between 0 and 180°!")

    # Validate arc
    if arc not in ['minor', 'major', 'both']:

      sys.exit ("Error: 'arc' must be 'minor' (clockwise)," +
                " 'major' (counter-clockwise), or 'both'!")

    # Validate min_loop
    if min_loop < 0:

      sys.exit ("Error: 'min_loop' must be a non-negative integer!")

    # Round source depth
    dp = int (round (source_depth))

    # Set file paths
    codes_path  = f"{self._tables_path}/codes_{dp}.list"
    phases_path = f"{self._tables_path}/phases_{dp}.list"
    table_path  = f"{self._tables_path}/table_{dp}.bin"

    # Initialize sets
    n_phases_set    = list ()
    phases_list_set = list ()

    Ni_set = list ()
    Rs_set = list ()
    Ds_set = list ()

    # Determine number of requested phases
    nr_phases = 0 if phases == 'all' else len (phases)

    # Prepare requested phases array
    requested_phases = np.zeros ((nr_phases, MAX_STRING_LEN), dtype = 'S1')    

    if nr_phases > 0:

      for i, p in enumerate (phases):
          
        # Fill phase names
        requested_phases[i, : len (p)] = list (p)

    for loop in range (min_loop, max_loop + 1):

      if arc == 'minor':

        # Compute minor delta
        delta_m = delta + float (loop) * 360.0

        n_phases_m = cp.c_int (0)
        
        phases_list_m = np.zeros ((MAX_PHASES, MAX_STRING_LEN), dtype = 'S1')

        ii_m = np.zeros (MAX_PHASES, dtype = np.int32)
        Ni_m = np.zeros (MAX_PHASES, dtype = np.int32)

        Rs_m = np.zeros ((MAX_PHASES, MAX_PATH_NODES), dtype = np.float64)
        Ds_m = np.zeros ((MAX_PHASES, MAX_PATH_NODES), dtype = np.float64)

        # Call ttimesAndPaths for minor
        return_value = libtracer.ttimesAndPaths (source_depth, delta_m,
                                                 cp.c_char_p (codes_path.encode ('utf-8')),
                                                 cp.c_char_p (phases_path.encode ('utf-8')),
                                                 cp.c_char_p (table_path.encode ('utf-8')),
                                                 self._ni.value,
                                                 self._nl.value,
                                                 self._rn, self._vpn, self._vsn,
                                                 self._nd, self._indices,
                                                 nr_phases,
                                                 requested_phases,
                                                 cp.byref (n_phases_m),
                                                 phases_list_m,
                                                 ii_m, Ni_m,
                                                 Rs_m, Ds_m)

        n_phases_loop = n_phases_m.value

        phases_list_loop = phases_list_m[:n_phases_loop]

        indices = ii_m[:n_phases_loop]
        Ni_loop = Ni_m[indices]
        
        Rs_loop = Rs_m[indices]
        Ds_loop = Ds_m[indices]

      elif arc == 'major':

        # Compute major delta
        delta_M = float (loop + 1) * 360.0 - delta

        n_phases_M = cp.c_int (0)
        
        phases_list_M = np.zeros ((MAX_PHASES, MAX_STRING_LEN), dtype = 'S1')

        ii_M = np.zeros (MAX_PHASES, dtype = np.int32)
        Ni_M = np.zeros (MAX_PHASES, dtype = np.int32)

        Rs_M = np.zeros ((MAX_PHASES, MAX_PATH_NODES), dtype = np.float64)
        Ds_M = np.zeros ((MAX_PHASES, MAX_PATH_NODES), dtype = np.float64)

        # Call ttimesAndPaths for major
        return_value = libtracer.ttimesAndPaths (source_depth, delta_M,
                                                 cp.c_char_p (codes_path.encode ('utf-8')),
                                                 cp.c_char_p (phases_path.encode ('utf-8')),
                                                 cp.c_char_p (table_path.encode ('utf-8')),
                                                 self._ni.value,
                                                 self._nl.value,
                                                 self._rn, self._vpn, self._vsn,
                                                 self._nd, self._indices,
                                                 nr_phases,
                                                 requested_phases,
                                                 cp.byref (n_phases_M),
                                                 phases_list_M,
                                                 ii_M, Ni_M,
                                                 Rs_M, Ds_M)

        n_phases_loop = n_phases_M.value
      
        phases_list_loop = phases_list_M[:n_phases_loop]

        indices = ii_M[:n_phases_loop]
        Ni_loop = Ni_M[indices]
        
        Rs_loop = Rs_M[indices]
        # Adjust deltas for major
        Ds_loop = np.pi - Ds_M[indices]

      elif arc == 'both':

        # Compute minor and major deltas
        delta_m = delta + float (loop) * 360.0
        delta_M = float (loop + 1) * 360.0 - delta

        # Min and max deltas
        delta_m, delta_M = min (delta_m, delta_M), max (delta_m, delta_M)

        n_phases_m = cp.c_int (0)
        n_phases_M = cp.c_int (0)

        phases_list_m = np.zeros ((MAX_PHASES, MAX_STRING_LEN), dtype = 'S1')
        phases_list_M = np.zeros ((MAX_PHASES, MAX_STRING_LEN), dtype = 'S1')

        ii_m = np.zeros (MAX_PHASES, dtype = np.int32)
        ii_M = np.zeros (MAX_PHASES, dtype = np.int32)

        Ni_m = np.zeros (MAX_PHASES, dtype = np.int32)
        Ni_M = np.zeros (MAX_PHASES, dtype = np.int32)

        Rs_m = np.zeros ((MAX_PHASES, MAX_PATH_NODES), dtype = np.float64)
        Rs_M = np.zeros ((MAX_PHASES, MAX_PATH_NODES), dtype = np.float64)

        Ds_m = np.zeros ((MAX_PHASES, MAX_PATH_NODES), dtype = np.float64)
        Ds_M = np.zeros ((MAX_PHASES, MAX_PATH_NODES), dtype = np.float64)

        # Call for minor
        return_value_m = libtracer.ttimesAndPaths (source_depth, delta_m,
                                                   cp.c_char_p (codes_path.encode ('utf-8')),
                                                   cp.c_char_p (phases_path.encode ('utf-8')),
                                                   cp.c_char_p (table_path.encode ('utf-8')),
                                                   self._ni.value,
                                                   self._nl.value,
                                                   self._rn, self._vpn, self._vsn,
                                                   self._nd, self._indices,
                                                   nr_phases,
                                                   requested_phases,
                                                   cp.byref (n_phases_m),
                                                   phases_list_m,
                                                   ii_m, Ni_m,
                                                   Rs_m, Ds_m)

        # Call for major
        return_value_M = libtracer.ttimesAndPaths (source_depth, delta_M,
                                                   cp.c_char_p (codes_path.encode ('utf-8')),
                                                   cp.c_char_p (phases_path.encode ('utf-8')),
                                                   cp.c_char_p (table_path.encode ('utf-8')),
                                                   self._ni.value,
                                                   self._nl.value,
                                                   self._rn, self._vpn, self._vsn,
                                                   self._nd, self._indices,
                                                   nr_phases,
                                                   requested_phases,
                                                   cp.byref (n_phases_M),
                                                   phases_list_M,
                                                   ii_M, Ni_M,
                                                   Rs_M, Ds_M)

        # Determine return value
        return_value = (return_value_m if return_value_m != 0
          else return_value_M if return_value_M != 0
          else 0
        )

        if return_value == 0:

          n_phases_m = n_phases_m.value
          n_phases_M = n_phases_M.value

          indices_m = ii_m[:n_phases_m]
          indices_M = ii_M[:n_phases_M]

          Ni_m = Ni_m[indices_m]
          Ni_M = Ni_M[indices_M]
        
          Rs_m = Rs_m[indices_m]
          Rs_M = Rs_M[indices_M]
          Ds_m = Ds_m[indices_m]
          # Adjust major deltas
          Ds_M = np.pi - Ds_M[indices_M]

          # Merge minor and major
          return_tuple = self._mergeMinorAndMajor (n_phases_m, phases_list_m,
                                                   Ni_m, Rs_m, Ds_m,
                                                   n_phases_M, phases_list_M,
                                                   Ni_M, Rs_M, Ds_M)

          n_phases_loop, phases_list_loop, Ni_loop, Rs_loop, Ds_loop = return_tuple

      if n_phases_loop == 0:
        
        # Warn if no phases
        warnings.warn (f"Warning: no phases found for loop {loop}!")

        continue

      if return_value == 1:

        sys.exit (f"Error: could not open file '{codes_path}'!")

      elif return_value == 2:

        sys.exit (f"Error: could not open file '{phases_path}'!")

      elif return_value == 3:

        sys.exit (f"Error: could not open file '{table_path}'!")

      elif return_value == 4:

        sys.exit (f"Error: could not read file '{codes_path}'!")

      elif return_value == 5:

        sys.exit (f"Error: could not read file '{phases_path}'!")

      elif return_value == 6:

        sys.exit (f"Error: could not read file '{table_path}'!")

      # Append to sets
      n_phases_set.append (n_phases_loop)

      phases_list_set.append (phases_list_loop)
      
      Ni_set.append (Ni_loop)
      
      Rs_set.append (Rs_loop)
      Ds_set.append (Ds_loop)

    if not n_phases_set:

      sys.exit ("Error: no phases found for the provided parameters!")

    # Initialize with first loop
    n_phases = n_phases_set[0]
    
    phases_list = phases_list_set[0]

    Ni = Ni_set[0]
    
    Rs = Rs_set[0]
    Ds = Ds_set[0]

    for loop in range (len (n_phases_set) - 1):

      # Merge subsequent loops
      n_phases, phases_list, Ni, Rs, Ds = self._mergePhases (n_phases,
                                                             phases_list,
                                                             Ni,
                                                             Rs,
                                                             Ds,
                                                             n_phases_set[loop + 1],
                                                             phases_list_set[loop + 1],
                                                             Ni_set[loop + 1],
                                                             Rs_set[loop + 1],
                                                             Ds_set[loop + 1])

    # Initialize phases dict
    phases_dict = dict ()

    if compute_amplitudes:

      # Initialize coefficient arrays
      cf_FP = np.empty ((self._N_TOAG, self._N_COEF), dtype = np.complex128)
      cf_FS = np.empty ((self._N_TOAG, self._N_COEF), dtype = np.complex128)
      cf_SP = np.empty ((self._N_TOAG, self._N_DISC, self._N_COES), dtype = np.complex128)
      cf_SS = np.empty ((self._N_TOAG, self._N_DISC, self._N_COES), dtype = np.complex128)

      # Compute source coefficients
      libcoefficients.computeCoefficients (source_depth, self._ni.value, self._nl.value,
                                           self._rn, self._rhon, self._vpn, self._vsn,
                                           self._nd, self._indices,
                                           cf_FP, cf_FS, cf_SP, cf_SS)

    for i in range (n_phases):

      # Parse info
      info = phases_list[i].tobytes ().split ()

      phase = info[0].decode ()

      if phase not in phases_dict:

        phases_dict[phase] = list ()

      # Parse values
      ttime    = float (info[1])
      take_off = float (info[2])
      dT_dD    = float (info[3])
      d2T_dD2  = float (info[4])
      
      code = info[5].decode ()

      # Get node count
      N = Ni[i]

      # Get input radii and deltas
      R_in = Rs[i][:N]
      D_in = Ds[i][:N]

      # Check boundary points
      is_boundary_point = self._isBoundaryPoint (R_in)

      # Compute extra nodes
      n_extra = int (is_boundary_point.sum ())

      # Initialize processed arrays
      R_proc = np.empty (N + n_extra, dtype = float)
      D_proc = np.empty (N + n_extra, dtype = float)

      j = 0

      for r, d, is_bp in zip (R_in, D_in, is_boundary_point):

        # Assign values
        R_proc[j] = r
        D_proc[j] = d
        
        j += 1

        if is_bp:
 
          # Duplicate for boundary
          R_proc[j] = r
          D_proc[j] = d
          
          j += 1

      # Get branches
      branches = self._branchesAlongPath (code, source_depth, R_proc, phase)

      # Compute extra counts
      extra_count_before = np.cumsum (is_boundary_point, dtype = np.int32) - \
                                      is_boundary_point.astype (np.int32)

      # Expand indices
      expanded_indices = np.arange (N, dtype = np.int32) + extra_count_before
      boundary_indices = np.flatnonzero (is_boundary_point)

      for index in boundary_indices:
        
        # Adjust boundary branches
        i0 = expanded_indices[index]; i1 = i0 + 1

        branches[i0] = branches[i0 - 1]
        branches[i1] = branches[i1 + 1] 

      # Create phase dict
      phase_dict = {'ttime' : ttime,
                    'take-off' : take_off,
                    'dT_dD' : dT_dD,
                    'd2T_dD2' : d2T_dD2,
                    'radii' : R_proc,
                    'deltas' : D_proc,
                    'branches' : branches,
                    'amplitude' : None}

      if compute_amplitudes:

        # Get radii and branches
        radii = phase_dict['radii']
        branches = phase_dict['branches']
        take_off = phase_dict['take-off']

        # Determine source type
        source_type = 'P'

        if branches[0] == 's' or branches[0] == 'S':
          
          source_type = 'S'

        # Compute toag index
        toag_index = int ((take_off / 90.0) * self._N_TOAG + 0.5)

        if toag_index < 0: toag_index = 0
        if toag_index > self._N_TOAG - 1: toag_index = self._N_TOAG - 1

        # Classify interactions
        interactions = self._classifyDiscontinuities (radii, branches)

        # Initialize cf
        cf = 1.0

        for interaction in interactions:

          disc = interaction[0]
          key  = interaction[1]

          disc_index = self._map[disc]
          
          if disc_index == 0:

            # Surface key index
            key_index = self._keysF[key]

            if source_type == 'P':

              cf *= cf_FP[toag_index, key_index]

            else:

              cf *= cf_FS[toag_index, key_index]

          else:

            # Internal key index
            key_index = self._keysS[key]

            if source_type == 'P':

              cf *= cf_SP[toag_index, disc_index - 1, key_index]

            else:

              cf *= cf_SS[toag_index, disc_index - 1, key_index]

        # Assign amplitude
        phase_dict['amplitude'] = cf

        # Append phase dict
        phases_dict[phase] += [phase_dict]

    # Create Phases object
    return Phases (source_depth, delta, phases_dict, self._rn[self._indices[:-1]])

  def _extendPath (self, ddelta, radius, radii, deltas, branches):
    """Extends a ray path at the turning point with a great circle at the CMB."""
    # Find min radius index
    imin = np.argmin (radii)
    # Compute average delta diff
    dd   = np.mean (np.abs (np.diff (deltas)))

    # Compute new sizes
    n       = len (deltas)
    n_extra = max (int (ddelta / dd + 0.5), 1)
    n_new   = n + n_extra

    # Initialize new arrays
    radii_new    = np.empty (n_new, dtype = np.float64)
    deltas_new   = np.empty (n_new, dtype = np.float64)
    branches_new = np.empty (n_new, dtype = '<U1')

    # Define indices
    i0 = 0
    i1 = imin
    i2 = imin + 1
    i3 = imin + n_extra
    i4 = n
    i5 = n_new

    # Assign radii sections
    radii_new[i0:i2] = radii[i0:i2]
    radii_new[i2:i3] = radius
    radii_new[i3]    = radii[i1]
    radii_new[i3:i5] = radii[i1:i4]

    # Get radius and delta values
    r0 = radii[imin - 1]
    r1 = radii[imin] 
    d0 = deltas[imin - 1]
    d1 = deltas[imin]
      
    # Compute slope
    dd_dr = (d1 - d0) / (r1 - r0)
    dr = r1 - radius
    dd_slope = dd_dr * dr

    # Compute min and max delta
    delta_min = deltas[i1] - dd_slope
    delta_max = deltas[i1] + dd_slope - ddelta

    # Assign deltas sections
    deltas_new[i0:i2] = deltas[i0:i2]
    deltas_new[i2:i3] = np.linspace (delta_min, delta_max, n_extra - 1)
    deltas_new[i3:i5] = deltas[i1:i4] - ddelta

    # Assign branches
    branches_new[:] = branches[0]

    return radii_new, deltas_new, branches_new   

  def _addDifractedWaves (self, source_depth = 0, delta = None,
                          phases = 'all', arc = 'minor'):
    """Adds difracted P and S waves to a list of phases and rays."""
    # Initialize requested phases
    requested_phases = list ()
    
    # Get CMB index
    index = self._indices[self._map['cmb']] - 1

    # Get radius and velocities
    radius = self._rn[index]
    alpha  = self._vpn[index]
    beta   = self._vsn[index]

    if phases == 'all':
      # Default phases
      requested_phases = ['P', 'S']

    elif 'Pdiff' in phases or 'Sdiff' in phases:

      if 'Pdiff' in phases:
        # Add P
        requested_phases.append ('P')

      if 'Sdiff' in phases:
        # Add S
        requested_phases.append ('S')

    # Initialize dict
    phases_dict = dict ()

    # Find deltas
    table = self._findDeltasForPhases (source_depth, requested_phases)

    if arc == 'minor':
      
      # Minor arc
      local_arcs = ['minor']

    elif arc == 'major':
      
      # Major arc
      local_arcs = ['major']

    else:

      # Both arcs
      local_arcs = ['minor', 'major']

    for local_arc in local_arcs:

      for phase_key in table:

        # Adjust max delta
        max_delta = table[phase_key][0, 0] - THRESHOLD
        # Compute phases
        phase = self._ttimesAndPaths (source_depth,
                                      max_delta,
                                      [phase_key])

        # Set phase name
        phase_name = 'Pdiff' if phase_key == 'P' else 'Sdiff'
        
        if arc == 'both':

          if local_arc == 'minor':

            phase_name += '_m'

          elif local_arc == 'major':

            phase_name += '_M'
        
        # Set velocity
        velocity   = alpha if phase_key == 'P' else beta

        # Compute local delta
        local_delta = delta if local_arc == 'minor' else 360.0 - delta
        ddelta      = np.radians (local_delta - max_delta)

        # Compute dT_dD and dttime
        dT_dD  = radius / velocity
        dttime = ddelta * dT_dD

        # Get phase dict
        phase_dict = phase.phases[phase_key][0]
        phase_dict['ttime'] += dttime

        # Get path components
        radii    = phase_dict['radii']
        deltas   = phase_dict['deltas']
        branches = phase_dict['branches']

        # Extend path
        radii, deltas, branches = self._extendPath (ddelta, radius,
                                                    radii, deltas,
                                                    branches)

        if local_arc == 'major':
          # Adjust deltas for major
          deltas = np.pi - deltas

        # Update dict
        phase_dict['deltas']    = deltas
        phase_dict['radii']     = radii
        phase_dict['branches']  = branches
        phase_dict['dT_dD']     = dT_dD
        phase_dict['d2T_dD2']   = 0.0
        phase_dict['amplitude'] = None

        # Add to phases dict
        phases_dict[phase_name] = [phase_dict]

    return phases_dict

  def ttimesAndPaths (self, source_depth = 0, delta = None, phases = 'all',
                      arc = 'minor', min_loop = 0, max_loop = 0,
                      compute_amplitudes = True):
    """Wrapper for _ttimesAndPaths and _addDifractedWaves."""
    if any (phase not in ['Pdiff', 'Sdiff'] for phase in phases):
      # Compute main phases
      phases_object = self._ttimesAndPaths (source_depth, delta,
                                            phases, arc,
                                            min_loop, max_loop,
                                            compute_amplitudes)
    
    else:
      # Initialize empty phases
      phases_object = Phases (source_depth, delta, dict (),
                              self._rn[self._indices[:-1]])

    if phases == 'all' or 'Pdiff' in phases or 'Sdiff' in phases:
      # Add diffracted waves
      phases_object.phases.update (self._addDifractedWaves (source_depth,
                                                            delta,
                                                            phases,
                                                            arc))

    return phases_object

  def getTravelTimes (self, source_depth = 0, phases = ['P', 'S']):
    """Get travel times of requested phases for a given source depth."""
    # Find deltas
    table = self._findDeltasForPhases (source_depth, phases)

    # Initialize ttimes dict
    ttimes = dict ()

    for phase in table:

      # Print computing message
      print (f"Computing travel times for {phase}...")

      if len (table[phase]) == 0:
        
        continue

      # Initialize phase dict
      ttimes[phase] = {'deltas' : [], 'ttimes' : [], 'dT_dD' : []}

      for delta, dT_dD in table[phase][:, :2]:

        try:

          # Ignore warnings
          warnings.simplefilter ('ignore')

          # Recover this row's arc and loop from its total distance
          d0   = delta % 360.0
          loop = int (delta // 360.0)

          if d0 <= 180.0:

            arc = 'minor'
            proper_delta = d0

          else:

            arc = 'major'
            proper_delta = 360.0 - d0

          # Compute arrivals
          arrivals = self.ttimesAndPaths (source_depth = source_depth,
                                          delta = proper_delta,
                                          phases = [phase], arc = arc,
                                          min_loop = loop, max_loop = loop)
          arrivals_dict = arrivals.getDict ()

          for arrival in arrivals_dict[phase]:

            # Convert ttime to minutes
            ttime = arrival['ttime'] / 60.0
            dT_dD = arrival['dT_dD']

            # Append values
            ttimes[phase]['deltas'].append (delta)
            ttimes[phase]['ttimes'].append (ttime)
            ttimes[phase]['dT_dD'].append (dT_dD)

        except:

          # Append nan for failure
          ttimes[phase]['deltas'].append (delta)
          ttimes[phase]['ttimes'].append (np.nan)
          ttimes[phase]['dT_dD'].append (dT_dD)

      if ttimes[phase]['deltas']:
  
        # Create arrays
        arrays = np.array ([ttimes[phase]['deltas'],
                            ttimes[phase]['ttimes'],
                            ttimes[phase]['dT_dD']
        ])
          
        # Sort by dT_dD
        sort_idx = np.argsort (arrays[2])
        sorted_arrays = arrays[:, sort_idx]

        y = sorted_arrays[1]
        p = sorted_arrays[2]

        # Valid mask
        valid_mask = ~np.isnan (y)

        if np.any (valid_mask):

          # Fill missing values by interpolating against dT_dD (single-valued)
          interpolator = interp1d (p[valid_mask], y[valid_mask],
                                   kind = 'linear',
                                   fill_value = 'extrapolate')

          y[np.isnan (y)] = interpolator (p)[np.isnan (y)]

        else:

          # Set to None if all nan
          y[:] = None

        # Break disjoint branches (negative-slope segments) so they are not
        # joined by a line; run after interpolation so breaks are not filled
        dl = sorted_arrays[0]
        tm = sorted_arrays[1]
        pr = sorted_arrays[2]

        step_min  = 0.05
        slope_tol = -0.01

        broken_d = [dl[0]]
        broken_t = [tm[0]]
        broken_p = [pr[0]]

        for k in range (1, dl.size):

          ddelta = dl[k] - dl[k - 1]
          dttime = tm[k] - tm[k - 1]

          if (abs (ddelta) >= step_min and np.isfinite (dttime) and
              dttime / ddelta < slope_tol):

            broken_d.append (np.nan)
            broken_t.append (np.nan)
            broken_p.append (np.nan)

          broken_d.append (dl[k])
          broken_t.append (tm[k])
          broken_p.append (pr[k])

        # Assign sorted lists (with branch breaks)
        ttimes[phase]['deltas'] = broken_d
        ttimes[phase]['ttimes'] = broken_t
        ttimes[phase]['dT_dD']  = broken_p

    return ttimes

  def plotTravelTimes (self, source_depth = 0, phases = ['P', 'S']):
    """Plot travel times of requested phases for a given source depth an."""
    # Get travel times
    ttimes = self.getTravelTimes (source_depth, phases)

    # Create figure
    plt.figure (figsize = (9, 7), dpi = 140)

    for phase in ttimes:

      # Plot curve
      plt.plot (ttimes[phase]['deltas'],
                ttimes[phase]['ttimes'],
                label = phase)

    # Set ticks
    plt.tick_params (axis = 'x', labelsize = 16)
    plt.tick_params (axis = 'y', labelsize = 16)
    # Set labels
    plt.ylabel ('Travel Time [minutes]', fontsize = 16)
    plt.xlabel (r'Delta [$\circ$]', fontsize = 16)
    # Set limits
    plt.xlim (0, None)
    plt.ylim (0, None)
    # Add grid
    plt.grid (linestyle = 'dashed', linewidth = '0.5')
    # Add legend
    plt.legend (fontsize = 12)
    # Set title
    plt.title ('Travel-Time Curves', fontsize = 20, y = 1.02)
    plt.tight_layout ()
    plt.show ()