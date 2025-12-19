/*
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

----------------------------------------------------------------------------------------------- */

#ifndef AMPLITUDES_H
#define AMPLITUDES_H
double complex freeSurface (enum keysF key,
                            double complex a, double complex b,
                            double complex c, double complex d);
double complex solid2Solid (enum keysS key,
                            double p, double psq,
                            double a, double b, double c, double d,
                            double alpha1, double alpha2,
                            double beta1, double beta2,
                            double rho1, double rho2,
                            double complex cosi1, double complex cosi2,
                            double complex cosj1, double complex cosj2,
                            double complex E, double complex F, double complex G,
                            double complex H, double complex D, double complex delta);
double complex solid2Liquid (enum keysS key,
                             double p, double a, double c,
                             double alpha1, double alpha2,
                             double beta1, double epsilon,
                             double complex cosi1, double complex cosi2,
                             double complex cosj1,
                             double complex D, double complex E, double complex F,
                             double complex delta);
double complex liquid2Solid (enum keysS key,
                             double p, double a, double c,
                             double alpha1, double alpha2,
                             double beta2, double epsilon,
                             double complex cosi1, double complex cosi2,
                             double complex cosj2,
                             double complex D, double complex E, double complex F,
                             double complex delta);
#endif