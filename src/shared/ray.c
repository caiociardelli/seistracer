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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "constants.h"
#include "exmath.h"
#include "enums.h"
#include "structs.h"

static double xpSph (double psph, double a, double b, double rb, double rt)
{
  /* Computes arc integral within a layer of semi-linear gradient. */
  /* Cast psph to complex */
  double complex p = psph;

  /* Compute ap and ap_rto */
  double complex ap = a * p;
  double complex ap_rto = csqrt (1.0 - csquare (ap) + EPSILON);

  /* Compute pv1 and pv2 */
  double complex pv1 = p * (a * rt + b);
  double complex pv2 = p * (a * rb + b);

  /* Compute k1_prdc and k2_prdc */
  double complex k1_prdc = -ap * pv1 + rt;
  double complex k2_prdc = -ap * pv2 + rb;

  /* Compute k1_sqrt and k2_sqrt */
  double complex k1_sqrt = csqrt (csquare (rt) - csquare (pv1) + EPSILON);
  double complex k2_sqrt = csqrt (csquare (rb) - csquare (pv2) + EPSILON);

  /* Compute k1_rto and k2_rto */
  double complex k1_rto = k1_prdc / (ap_rto * k1_sqrt + EPSILON);
  double complex k2_rto = k2_prdc / (ap_rto * k2_sqrt + EPSILON);

  /* Compute k1 and k2 */
  double complex k1 = (ap * catanh (k1_rto)) / (ap_rto + EPSILON)
                       - catan (pv1 / (k1_sqrt + EPSILON));
  double complex k2 = (ap * catanh (k2_rto)) / (ap_rto + EPSILON)
                       - catan (pv2 / (k2_sqrt + EPSILON));

  /* Return real part with a negative sign */
  return -creal (k2 - k1);
}

static double tpSph (double psph, double a, double b, double rb, double rt)
{
  /* Computes time integral within a layer of semi-linear gradient. */
  /* Handle small slope case */
  if (fabs (a) < MIN_SLOPE)
  {
    double complex bpsph = b * psph;

    double complex k1 = csqrt (square (rt) - csquare (bpsph))
                      / (b + EPSILON);
    double complex k2 = csqrt (square (rb) - csquare (bpsph))
                      / (b + EPSILON);

    return -creal (k2 - k1);
  }

  /* Cast psph to complex */
  double complex p = psph;

  /* Compute ap and ap_rto */
  double complex ap = a * p;
  double complex ap_rto = csqrt (1.0 - csquare (ap) + EPSILON);

  /* Compute pv1 and pv2 */
  double complex pv1 = p * (a * rt + b);
  double complex pv2 = p * (a * rb + b);

  /* Compute k1_prdc and k2_prdc */
  double complex k1_prdc = -ap * pv1 + rt;
  double complex k2_prdc = -ap * pv2 + rb;

  /* Compute k1_sqrt and k2_sqrt */
  double complex k1_sqrt = csqrt (csquare (rt) - csquare (pv1) + EPSILON);
  double complex k2_sqrt = csqrt (csquare (rb) - csquare (pv2) + EPSILON);

  /* Compute k1_rto and k2_rto */
  double complex k1_rto = k1_prdc / (ap_rto * k1_sqrt + EPSILON);
  double complex k2_rto = k2_prdc / (ap_rto * k2_sqrt + EPSILON);

  /* Compute k1 and k2 */
  double complex k1 = (catanh (k1_rto) / (ap_rto + EPSILON)
                       - catanh (rt / (k1_sqrt + EPSILON))) / (a + EPSILON);
  double complex k2 = (catanh (k2_rto) / (ap_rto + EPSILON)
                       - catanh (rb / (k2_sqrt + EPSILON))) / (a + EPSILON);

  /* Return real part with a negative sign */
  return -creal (k2 - k1);
}

static void downwardFast (int phase,
                          int l, int n,
                          int *m, int *L,
                          double psph, double rs,
                          double *d,
                          int nl, double rn[nl],
                          int nd, int Id[nd],
                          struct Model Md[n])
{
  /* Computes delta of a downgoing ray within a layer of semi-linear
     gradient. */
  /* Get layer boundaries */
  double zrb = rn[Id[l]], zrt = rn[Id[l - 1]];

  while (*m >= 0 && *m < n)
  {
    double rb = Md[*m].rb, rt = Md[*m].rt;

    if (rt > zrt || rb < zrb) break;

    /* Get velocities and coefficients */
    double vb, vt, a, b;

    switch (phase)
    {
      case P:
        vb = Md[*m].vpb; a = Md[*m].avp;
        vt = Md[*m].vpt; b = Md[*m].bvp;
      break;

      case S:
        vb = Md[*m].vsb; a = Md[*m].avs;
        vt = Md[*m].vst; b = Md[*m].bvs;
      break;

      default:
        fprintf (stderr, "Unknown phase.\n");
      break;
    }
  
    if (rt / vt < psph) break;
    if (*L == 1 && rt > rs) {rt = rs; vt = a * rs + b;}
    if (rb / vb < psph) rb = b / (1.0 / psph - a);
    
    rb = fmax (rb, Md[*m].rb);

    /* Update delta */
    *d -= xpSph (psph, a, b, rb, rt);

    *L += 1; *m += 1;

    if (*L > MAX_PATH_NODES)
    {
      fprintf (stderr, "Error: path has more nodes than MAX_PATH_NODES!\n");

      exit (EXIT_FAILURE);
    }
  }
}

static void upwardFast (int phase,
                        int l, int n,
                        int *m, int *L,
                        double psph, double rs,
                        double *d,
                        int nl, double rn[nl],
                        int nd, int Id[nd],
                        struct Model Md[n])
{
  /* Computes delta of a upgoing ray within a layer of semi-linear
     gradient. */
  /* Get layer boundaries */
  double zrb = rn[Id[l]], zrt = rn[Id[l - 1]];

  while (*m >= 0)
  {
    double rb = Md[*m].rb, rt = Md[*m].rt;

    if (rt > zrt || rb < zrb) break;

    /* Get velocities and coefficients */
    double vb, a, b;

    switch (phase)
    {
      case P:
        vb = Md[*m].vpb; a = Md[*m].avp;
                         b = Md[*m].bvp;
      break;

      case S:
        vb = Md[*m].vsb; a = Md[*m].avs;
                         b = Md[*m].bvs;
      break;

      default:
        fprintf (stderr, "Unknown phase.\n");
      break;
    }

    if (*L == 1 && rs > rb) {rb = rs; vb = a * rs + b;}
    if (rb / vb < psph) rb = b / (1.0 / psph - a);
    
    rb = fmax (rb, Md[*m].rb);

    /* Update delta */
    *d -= xpSph (psph, a, b, rb, rt);

    *L += 1; *m -= 1;

    if (*L > MAX_PATH_NODES)
    {
      fprintf (stderr, "Error: path has more nodes than MAX_PATH_NODES!\n");

      exit (EXIT_FAILURE);
    }
  }
}

double computeDelta (int l, int m,
                     int n, struct Model Md[n],
                     int np, int Pm[MAX_HEAP_DEPTH],
                     int nl, double rn[nl],
                     int nd, int Id[nd],
                     double psph,
                     double rs, double ds)
{
  /* Computes delta. */
  if (psph < MIN_PSPH) psph = MIN_PSPH;

  /* Initialize L and d */
  int L = 1; double d = ds;

  for (int i = 0; i < np; i++)
  {
    if (i > 0)
    {
      /* Adjust layer and model index */
      if ((Pm[i - 1] == Pd ||
           Pm[i - 1] == Sd) &&
          (Pm[i] == Pu ||
           Pm[i] == Su)) {l--; m--;}

      else if ((Pm[i - 1] == Pu ||
                Pm[i - 1] == Su) &&
               (Pm[i] == Pd ||
                Pm[i] == Sd)) {l++; m++;}
    }

    switch (Pm[i])
    {
      case Pu:
        upwardFast (P, l, n, &m, &L, psph, rs, &d, nl, rn, nd, Id, Md); l--;
      break;

      case Su:
        upwardFast (S, l, n, &m, &L, psph, rs, &d, nl, rn, nd, Id, Md); l--;
      break;

      case Pd:
        downwardFast (P, l, n, &m, &L, psph, rs, &d, nl, rn, nd, Id, Md); l++;
      break;

      case Sd:
        downwardFast (S, l, n, &m, &L, psph, rs, &d, nl, rn, nd, Id, Md); l++;
      break;

      default:
        fprintf (stderr, "Unknown phase.\n");
      break;
    }
  }

  /* Return degree conversion */
  return rad2Degree (fabs (d - ds));
}

static void downward (int phase,
                      int l, int n,
                      int *m, int *L,
                      double psph, double rs,
                      double *d, double *t,
                      double Rs[MAX_PATH_NODES],
                      double Ds[MAX_PATH_NODES],
                      int nl, double rn[nl],
                      int nd, int Id[nd],
                      struct Model Md[n])
{
  /* Computes path of downgoing ray within a layer of semi-linear
     gradient. */
  /* Get layer boundaries */
  double zrb = rn[Id[l]], zrt = rn[Id[l - 1]];

  while (*m >= 0 && *m < n)
  {
    double rb = Md[*m].rb, rt = Md[*m].rt;

    if (rt > zrt || rb < zrb) break;

    /* Get velocities and coefficients */
    double vb, vt, a, b;

    switch (phase)
    {
      case P:
        vb = Md[*m].vpb; a = Md[*m].avp;
        vt = Md[*m].vpt; b = Md[*m].bvp;
      break;

      case S:
        vb = Md[*m].vsb; a = Md[*m].avs;
        vt = Md[*m].vst; b = Md[*m].bvs;
      break;

      default:
        fprintf (stderr, "Unknown phase.\n");
      break;
    }

    if (rt / vt < psph) break;
    if (*L == 1 && rt > rs) {rt = rs; vt = a * rs + b;}
    if (rb / vb < psph) rb = b / (1.0 / psph - a);
    
    rb = fmax (rb, Md[*m].rb);

    /* Update d and t */
    *d -= xpSph (psph, a, b, rb, rt);
    *t += tpSph (psph, a, b, rb, rt);

    /* Store path data */
    Rs[*L] = rb; Ds[*L] = *d; *L += 1; *m += 1;

    if (*L > MAX_PATH_NODES)
    {
      fprintf (stderr, "Error: path has more nodes than MAX_PATH_NODES!\n");

      exit (EXIT_FAILURE);
    }
  }
}

static void upward (int phase,
                    int l, int n,
                    int *m, int *L,
                    double psph, double rs,
                    double *d, double *t,
                    double Rs[MAX_PATH_NODES],
                    double Ds[MAX_PATH_NODES],
                    int nl, double rn[nl],
                    int nd, int Id[nd],
                    struct Model Md[n])
{
  /* Computes path of upgoing ray within a layer of semi-linear
     gradient. */
  /* Get layer boundaries */
  double zrb = rn[Id[l]], zrt = rn[Id[l - 1]];

  while (*m >= 0)
  {
    double rb = Md[*m].rb, rt = Md[*m].rt;

    if (rt > zrt || rb < zrb) break;

    /* Get velocities and coefficients */
    double vb, a, b;

    switch (phase)
    {
      case P:
        vb = Md[*m].vpb; a = Md[*m].avp;
                         b = Md[*m].bvp;
      break;

      case S:
        vb = Md[*m].vsb; a = Md[*m].avs;
                         b = Md[*m].bvs;
      break;

      default:
        fprintf (stderr, "Unknown phase.\n");
      break;
    }

    if (*L == 1 && rs > rb) {rb = rs; vb = a * rs + b;}
    if (rb / vb < psph) rb = b / (1.0 / psph - a);
    
    rb = fmax (rb, Md[*m].rb);

    /* Update d and t */
    *d -= xpSph (psph, a, b, rb, rt);
    *t += tpSph (psph, a, b, rb, rt);

    /* Store path data */
    Rs[*L] = rt; Ds[*L] = *d; *L += 1; *m -= 1;

    if (*L > MAX_PATH_NODES)
    {
      fprintf (stderr, "Error: path has more nodes than MAX_PATH_NODES!\n");

      exit (EXIT_FAILURE);
    }
  }
}

void computeDeltaAndTime (int l, int m, int *L,
                          int n, struct Model Md[n],
                          int np, int Pm[MAX_HEAP_DEPTH],
                          int nl, double rn[nl],
                          int nd, int Id[nd],
                          double psph, double *delta, double *ttime,
                          double rs, double ds, double ts,
                          double Rs[MAX_PATH_NODES],
                          double Ds[MAX_PATH_NODES])
{
  /* Computes delta, ray path and travel time. */
  if (psph < MIN_PSPH) psph = MIN_PSPH;

  *L = 1; Rs[0] = rs; Ds[0] = ds;

  double d = ds, t = ts;

  for (int i = 0; i < np; i++)
  {
    if (i > 0)
    {
      /* Adjust layer for direction change */
      if ((Pm[i - 1] == Pd ||
           Pm[i - 1] == Sd) &&
          (Pm[i] == Pu ||
           Pm[i] == Su)) {l--; m--;}

      else if ((Pm[i - 1] == Pu ||
                Pm[i - 1] == Su) &&
               (Pm[i] == Pd ||
                Pm[i] == Sd)) {l++; m++;}
    }

    switch (Pm[i])
    {
      case Pu:
        upward (P, l, n, &m, L, psph, rs,
                &d, &t, Rs, Ds, nl, rn,
                nd, Id, Md); l--;
      break;

      case Su:
        upward (S, l, n, &m, L, psph, rs,
                &d, &t, Rs, Ds, nl, rn,
                nd, Id, Md); l--;
      break;

      case Pd:
        downward (P, l, n, &m, L, psph, rs,
                  &d, &t, Rs, Ds, nl, rn,
                  nd, Id, Md); l++;
      break;

      case Sd:
        downward (S, l, n, &m, L, psph, rs,
                  &d, &t, Rs, Ds, nl, rn,
                  nd, Id, Md); l++;
      break;

      default:
        fprintf (stderr, "Unknown phase.\n");
      break;
    }
  }

  *delta = rad2Degree (fabs (d - ds));
  *ttime = t - ts;
}