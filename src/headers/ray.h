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

#ifndef RAY_H
#define RAY_H
double computeDelta (int l, int m,
                     int n, struct Model Md[n],
                     int np, int Pm[MAX_HEAP_DEPTH],
                     int nl, double rn[nl],
                     int nd, int Id[nd],
                     double psph,
                     double rs, double ds);
void computeDeltaAndTime (int l, int m, int *L,
                          int n, struct Model Md[n],
                          int np, int Pm[MAX_HEAP_DEPTH],
                          int nl, double rn[nl],
                          int nd, int Id[nd],
                          double psph, double *delta, double *ttime,
                          double rs, double ds, double ts,
                          double Rs[MAX_PATH_NODES],
                          double Ds[MAX_PATH_NODES]);
#endif