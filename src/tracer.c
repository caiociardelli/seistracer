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
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include "constants.h"
#include "io.h"
#include "exmath.h"
#include "enums.h"
#include "structs.h"
#include "enums.h"
#include "source_and_model.h"
#include "heapsort.h"
#include "ray.h"
#include <string.h>

static void adjustDelta (double *delta)
{
  /* Adjusts the delta to avoid exact 180-degree multiples. */
  double k = round (*delta / 180.0);
  double target = k * 180.0;
  double diff = *delta - target;

  if (fabs (diff) < MAX_ERROR)
  {
    if      (diff < 0.0) *delta = target - MAX_ERROR;
    else if (diff > 0.0) *delta = target + MAX_ERROR;
    else
    {
      /* Handle exact multiple based on parity of k. */
      int ik = (int) k;

      if (ik % 2 == 0)
      
        *delta = target + MAX_ERROR;
      
      else
      
        *delta = target - MAX_ERROR;
    }
  }
}

static void codeToPermutation (char code[MAX_STRING_LEN], int N, int P[N])
{
  /* Computes the permutation array from the corresponding ascii code. */
  for (int i = 0; i < N; i ++)

    switch (code[i])
    {
      case 'A':
        P[i] = Pu;
      break;

      case 'B':
        P[i] = Su;
      break;

      case 'C':
        P[i] = Pd;
      break;

      case 'D':
        P[i] = Sd;
      break;

      default:
        fprintf (stderr, "Unknown permutation key %c.\n", code[i - 1]);
      break;
    }
}

static bool isNotInList (int index, int n, int list[n])
{
  /* Checks whether an index is on the list. */
  for (int i = 0; i < n; i ++)

    if (list[i] == index) return false;

  return true;
}

int ttimesAndPaths (double depth, double delta,
                    char codes_path[MAX_STRING_LEN],
                    char phases_path[MAX_STRING_LEN],
                    char table_path[MAX_STRING_LEN],
                    int ni,
                    int nl,
                    double rn[nl],
                    double vpn[nl],
                    double vsn[nl],
                    int nd, int Id[nd],
                    int nr_phases,
                    char requested_phases[nr_phases][MAX_STRING_LEN],
                    int *n_phases,
                    char phases_list[MAX_PHASES][MAX_STRING_LEN],
                    int ii[MAX_PHASES], int Ni[MAX_PHASES],
                    double Rs[MAX_PHASES][MAX_PATH_NODES],
                    double Ds[MAX_PHASES][MAX_PATH_NODES])
{
  /* Computes travel-times and paths of all possibles phases
    given the source depth and delta. */
  double odelta = delta;

  /* Adjust delta value */
  adjustDelta (&delta);
  
  /* Compute source radius and velocities */
  double rs = depth2R (depth), ds = 0.5 * PI, ts = 0.0;
  double vsP, vsS;

  sourceVelocity (rs, nl, rn, vpn, vsn, &vsP, &vsS);

  /* Create model structure array */
  struct Model M[ni];

  arraysToStructure (ni, M, nl, rn, vpn, vsn);

  /* Find layer indices */
  int li = setLayer (rs, nl, rn, nd, Id), mi = 0;

  while (M[mi].rb >= rs) mi++;

  /* Open input files */
  FILE *file1 = fopen (codes_path, "r");
  FILE *file2 = fopen (phases_path, "r");
  FILE *file3 = fopen (table_path, "r");

  if (file1 == NULL) return 1;
  if (file2 == NULL) return 2;
  if (file3 == NULL) return 3;

  /* Read number of codes */
  int nC, max_code_len = MAX_HEAP_DEPTH + 1;

  if (fscanf (file1, "%d", &nC) != 1) return 4;

  /* Allocate codes array */
  char (*codes)[max_code_len] = malloc (sizeof (char[nC][max_code_len]));

  if (codes == NULL)
  {
    fprintf (stderr, "Error: could not allocate memory for"
                     " codes...\n"); exit (EXIT_FAILURE);
  }

  /* Read codes */
  for (int i = 0; i < nC; i++)

    if (fscanf (file1, "%s", codes[i]) == 0) return 4;

  /* Read number of phases */
  int nP; if (fscanf (file2, "%d", &nP) != 1) return 5;

  /* Allocate Nl and phases */
  int Nl[nP];

  char (*phases)[MAX_STRING_LEN] = malloc (sizeof (char[nP][MAX_STRING_LEN]));

  if (phases == NULL)
  {
    fprintf (stderr, "Error: could not allocate memory for"
                     " phases...\n"); exit (EXIT_FAILURE);
  }

  /* Read phase data */
  int indices[nP], nvi = 0;

  for (int i = 0; i < nP; i++)
  {
    if (fscanf (file2, "%d %s", &Nl[i], phases[i]) != 2) return 5;

    for (int j = 0; j < nr_phases; j++)

      if (strcmp (phases[i], requested_phases[j]) == 0)

        indices[nvi++] = i;
  }

  /* Allocate string array */
  char (*string)[MAX_STRING_LEN] = malloc (sizeof (char[MAX_PHASES][MAX_STRING_LEN]));

  if (string == NULL)
  {
    fprintf (stderr, "Error: could not allocate memory for"
                     " string...\n"); exit (EXIT_FAILURE);
  }

  /* Read number of blocks */
  int N; if (fread (&N, sizeof (int), 1, file3) != 1) return 6;

  /* Allocate blocks */
  struct Block *blocks = malloc (sizeof (struct Block[N]));

  if (blocks == NULL)
  {
    fprintf (stderr, "Error: could not allocate memory for"
                     " blocks...\n"); exit (EXIT_FAILURE);
  }

  /* Read blocks data */
  if (fread (blocks, sizeof (struct Block), N, file3) != (size_t) N) return 6;

  /* Initialize counters and variables */
  int n = 0, c = 0, L, np1, np2;

  double psph1, psph2;
  double delta1, delta2;
  double tt[MAX_PHASES];

  for (int i = 0; i < nP; i++)
  {
    if (nr_phases > 0 && isNotInList (i, nvi, indices))
    {
      n += Nl[i]; continue;
    }

    for (int j = 0; j < Nl[i]; j++)
    {
      /* Read block data */
      int np = (int) blocks[n].np;
      int index = blocks[n].index;

      double delta_r = (double) blocks[n].delta;
      double psph_r  = (double) blocks[n++].psph;

      if (j == 0)
      {
        delta1 = delta_r; psph1 = psph_r; np1 = np;

        continue;
      }

      np2 = np; psph2 = psph_r; delta2 = delta_r;

      if (((delta > delta1 && delta <= delta2) ||
           (delta < delta1 && delta >= delta2)) && np1 == np2)
      {
        /* Create permutation array */
        int Pm[np];

        codeToPermutation (codes[index], np, Pm);

        double delta_i = delta, ttime, psph, error;

        double psph_i1 = psph1;
        double psph_i2 = psph2;

        double delta_i1 = delta1;
        double delta_i2 = delta2;

        /* Iterate to find psph */
        for (int k = 0; k < MAX_ITERATIONS; k++)
        {
          double a = (psph_i2 - psph_i1)
                   / (delta_i2 - delta_i1);
          double b = psph_i1 - a * delta_i1;

          psph = a * delta + b;

          delta_i = computeDelta (li, mi, ni, M, np, Pm, nl, rn, nd, Id, psph, rs, ds);

          if ((delta_i2 - delta_i1) * (delta_i - delta) < 0.0)
          {
            psph_i1 = psph; delta_i1 = delta_i;
          }

          else
          {
            psph_i2 = psph; delta_i2 = delta_i;
          }

          error = fabs (delta_i - delta);

          if (error < MIN_ERROR) break;
        }

        if (error < MAX_ERROR)
        {
          /* Determine source velocity */
          double vs = (Pm[0] == Pd || Pm[0] == Pu) ? vsP : vsS;
          /* Compute take-off angle */
          double take_off = rad2Degree (asin (vs * psph / rs));
          /* Compute psph derivative */
          double psph_prime = (psph_i2 - psph_i1) / (delta_i2 - delta_i1);

          /* Compute delta and time */
          computeDeltaAndTime (li, mi, &L, ni, M, np, Pm, nl, rn, nd, Id, psph,
                               &delta_i, &ttime, rs, ds, ts, Rs[c], Ds[c]);

          /* Rescale deltas */
          double rescaling = odelta / delta_i;

          for (int l = 0; l < L; l++)

            Ds[c][l] = (Ds[c][l] - ds) * rescaling + ds;

          /* Format output string */
          sprintf (string[c], "%-45s %11.2lf %11.3lf %12.3lf %10.3lf %s\n", phases[i],
                                                                            round (ttime * 100.0) / 100.0,
                                                                            take_off,
                                                                            psph,
                                                                            psph_prime,
                                                                            codes[index]);

          ii[c] = c; tt[c] = ttime; Ni[c++] = L;

          if (c == MAX_PHASES)
          {
            fprintf (stderr, "Error: Number of phases larger than MAX_PHASES!\n");

            exit (EXIT_FAILURE);
          }
        }
      }

      np1 = np2; delta1 = delta2; psph1 = psph2;
    }
  }

  /* Close files */
  fclose (file1);
  fclose (file2);
  fclose (file3);

  /* Sort by travel time */
  heapSort (c, ii, tt);

  *n_phases = c;

  for (int i = 0; i < c; i++)

    sprintf (phases_list[i], "%s", string[ii[i]]);

  /* Free allocated memory */
  free (codes);
  free (phases);
  free (blocks);
  free (string);

  return 0;
}