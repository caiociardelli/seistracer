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

#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H
void computeCoefficients (double depth,
                          int ni,
                          int nl,
                          double rn[nl],
                          double rhon[nl],
                          double vpn[nl],
                          double vsn[nl],
                          int nd, int Id[nd],
                          double complex cf_FP[N_TOAG][N_COEF],
                          double complex cf_FS[N_TOAG][N_COEF],
                          double complex cf_SP[N_TOAG][N_DISC][N_COES],
                          double complex cf_SS[N_TOAG][N_DISC][N_COES])
#endif