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

#ifndef SOURCE_AND_MODEL_H
#define SOURCE_AND_MODEL_H
void sourceVelocity (double rs,
                     int nl,
                     double rn[nl],
                     double vpn[nl],
                     double vsn[nl],
                     double *vs_P, double *vs_S);

int setLayer (double rs,
              int nl, double rn[nl],
              int nd, int Id[nd]);

void setInterfaces (int nl,
                    double rn[nl],
                    double rhon[nl],
                    double vpn[nl],
                    double vsn[nl],
                    int nd, int Id[nd],
                    double Psph_BP[N_DISC], double Psph_BS[N_DISC],
                    struct Surface *sfc,
                    struct Discontinuity dsc[N_DISC]);
void arraysToStructure (int n,
                        struct Model Md[n],
                        int nl,
                        double rn[nl],
                        double vpn[nl],
                        double vsn[nl]);
#endif