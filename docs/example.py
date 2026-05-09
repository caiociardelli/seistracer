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

import seistracer


if __name__ == '__main__':

  tracer = seistracer.Tracer (model = 'iasp91')

  print ("\nPlotting velocity model...")

  tracer.plotModel ()

  print ("\nComputing phases 'P' and 'S'...")

  phases = tracer.ttimesAndPaths (source_depth = 0, delta = 20, phases = ['P', 'S'])

  print ('\nInfo for all phases:\n')

  phases.getInfo ()

  print ('\nPlotting path for all phases...')

  phases.plotPath ()

  print ('\nPlotting travel-time curves (default are P and S)...')

  tracer.plotTravelTimes ()

  print ("\nComputing phases 'Pdiff' and 'Sdiff'...")

  phases = tracer.ttimesAndPaths (source_depth = 0, delta = 150, phases = ['Pdiff', 'Sdiff'])

  print ('\nInfo for all phases:\n')

  phases.getInfo ()

  print ('\nPlotting path for all phases...')

  phases.plotPath ()

  print ('\nPlotting travel-time curves for P, S, Pdiff, and Sdiff...')

  tracer.plotTravelTimes (phases = ['P', 'S', 'Pdiff', 'Sdiff'])

  print ('\nPlotting travel-time curves for multiple phases...')

  tracer.plotTravelTimes (phases = ['PP', 'SS', 'Pdiff', 'Sdiff', 'PS', 'SKP', 
                                    'SKKS', 'PKIKP', 'PKJKP', 'SKS660t+P'])

  print ("\nPlotting amplitude coefficients for all discontinuities...")

  tracer.plotAmplitudeCoefficients (discontinuities = 'all')

  print ("\nComputing phases 'P', 'P410-P', and 'SKS660t+P'...")

  phases = tracer.ttimesAndPaths (source_depth = 0, delta = 98,
                                  phases = ['P', 'P410-P', 'SKS660t+P'])

  print ('\nPlotting path for all three phases...')

  phases.plotPath ()

  print ('\nInfo for all three phases:\n')

  phases.getInfo ()

  print ("\nComputing all available phases...")

  phases = tracer.ttimesAndPaths (source_depth = 480, delta = 140, arc = 'minor')

  print ('\nPlotting path for all minor-arc phases...')

  phases.plotPath ()

  phases = tracer.ttimesAndPaths (source_depth = 480, delta = 140, arc = 'major')

  print ('\nPlotting path for all major-arc phases...')

  phases.plotPath ()

  phases = tracer.ttimesAndPaths (source_depth = 480, delta = 140, arc = 'both',
                                  compute_amplitudes = True)

  print ('\nPlotting path for all minor-arc and major-arc phases...')

  phases.plotPath ()

  print ('\nInfo for all phases:\n')

  phases.getInfo ()

  print ('\nInfo for some phases:\n')

  phases.getInfo (['pPKJKS', 'P410-ScSScS'])

  print ('\nPrint path for a single phase:\n')

  paths = phases.getPath (['pPKJKS'])
  print (paths)

  print ('\nGet all the information for a phase as a dictionary:\n')

  phases_dict = phases.getDict (['pPKJKS'])
  print (phases_dict)

  print ('\nPlotting path for some phases...')

  phases.plotPath (['PKIIKP', 'PKP660t+SScP', 'PKJKS',
                    'SPS660-S', 'ScSScS660-ScS_M', 'pPKJKS',
                    'pPmt-ScSScS', 'P410-ScSScS', 'sSKIKPPm+P'])

  print ('\nPlotting path for all minor-arc phases including the higher-order ones...')

  phases = tracer.ttimesAndPaths (source_depth = 480, delta = 140, arc = 'minor',
                                  min_loop = 0, max_loop = 2)

  phases.plotPath ()

  print ("\nComputing all available phases...")
  
  phases = tracer.ttimesAndPaths (source_depth = 20, delta = 150)

  print ('\nInfo for all phases:\n')

  phases.getInfo ()