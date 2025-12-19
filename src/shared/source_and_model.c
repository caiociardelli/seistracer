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

#include <math.h>
#include "constants.h"
#include "structs.h"

void sourceVelocity (double rs,
                     int nl,
                     double rn[nl],
                     double vpn[nl],
                     double vsn[nl],
                     double *vs_P, double *vs_S)
{
  /* Computes P- and S-wave velocities at the source. */
  /* Find layer index */
  int i = 1;

  while (rn[i] >= rs)
  {
    i++; if (i == nl) break;
  }

  /* Get layer boundaries and velocities */
  double rt  =  rn[i - 1], rb  =  rn[i];
  double vpt = vpn[i - 1], vpb = vpn[i];
  double vst = vsn[i - 1], vsb = vsn[i];

  double a, b;

  /* Compute P velocity coefficients */
  a = (vpt - vpb) / (rt - rb);
  b = vpb - a * rb;

  *vs_P = a * rs + b;

  /* Compute S velocity coefficients */
  a = (vst - vsb) / (rt - rb);
  b = vsb - a * rb;

  *vs_S = a * rs + b;
}

int setLayer (double rs,
              int nl, double rn[nl],
              int nd, int Id[nd])
{
  /* Sets starting index according to the layer
    in which the hypocenter is. */
  /* Check surface or crust */
  if (rs == rn[Id[0]] || (rs >= rn[Id[1]] &&
                          rs <  rn[Id[0]])) return 1;
  /* Check upper mantle */
  else if (rs >= rn[Id[2]] && rs < rn[Id[1]]) return 2;
  /* Check transition zone */
  else if (rs >= rn[Id[3]] && rs < rn[Id[2]]) return 3;
  /* Default to lower mantle */
  else return 4;
}

void setInterfaces (int nl,
                    double rn[nl],
                    double rhon[nl],
                    double vpn[nl],
                    double vsn[nl],
                    int nd, int Id[nd],
                    double Psph_BP[N_DISC], double Psph_BS[N_DISC],
                    struct Surface *sfc,
                    struct Discontinuity dsc[N_DISC])
{
  /* Extracts model parameters at surface and internal discontinuities. */
  /* Set surface properties */
  sfc->r     = rn[Id[0]];
  sfc->rho   = rhon[Id[0]];
  sfc->alpha = vpn[Id[0]];
  sfc->beta  = vsn[Id[0]];

  for (int l = 1; l <= N_DISC; l++)
  {
    /* Compute bottom psph for P and S */
    Psph_BP[l - 1] = rn[Id[l]] / vpn[Id[l] - 1];
    Psph_BS[l - 1] = rn[Id[l]] / vsn[Id[l] - 1];

    /* Set discontinuity properties */
    dsc[l - 1].r        = rn[Id[l]];
    dsc[l - 1].rho1     = rhon[Id[l] - 1];
    dsc[l - 1].rho2     = rhon[Id[l]];
    dsc[l - 1].alpha1   = vpn[Id[l] - 1];
    dsc[l - 1].alpha2   = vpn[Id[l]];
    dsc[l - 1].beta1    = vsn[Id[l] - 1];
    dsc[l - 1].beta2    = vsn[Id[l]];
  }
}

void arraysToStructure (int n,
                        struct Model Md[n],
                        int nl,
                        double rn[nl],
                        double vpn[nl],
                        double vsn[nl])
{
  /* Converts arrays to structure of arrays. */
  for (int l = 1, m = 0; l < nl; l++)

    if (fabs (rn[l] - rn[l - 1]) > EPSILON)
    {
      /* Set radii */
      Md[m].rb  =  rn[l]; Md[m].rt  =  rn[l - 1];
      /* Set velocities */
      Md[m].vpb = vpn[l]; Md[m].vpt = vpn[l - 1];
      Md[m].vsb = vsn[l]; Md[m].vst = vsn[l - 1];

      /* Compute P coefficients */
      double a, b, dr_1 = 1.0 / (Md[m].rt - Md[m].rb);

      a = (Md[m].vpt - Md[m].vpb) * dr_1;
      b = Md[m].vpb - a * Md[m].rb;

      Md[m].avp = a; Md[m].bvp = b;

      /* Compute S coefficients */
      a = (Md[m].vst - Md[m].vsb) * dr_1;
      b = Md[m].vsb - a * Md[m].rb;

      Md[m].avs = a; Md[m].bvs = b;

      m++;
    }
}