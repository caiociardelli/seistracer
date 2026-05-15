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
#include <locale.h>
#include "constants.h"
#include "constants_export.h"

void getSeisTracerConstants (struct SeisTracerConstants *out)
{
  /* Exports selected build-time constants to Python (ctypes). */
  if (out == NULL) return;

  /* Assign constant values */
  out->max_string_len = MAX_STRING_LEN;
  out->max_phases     = MAX_PHASES;
  out->max_path_nodes = MAX_PATH_NODES;
  out->n_toag         = N_TOAG;
  out->n_coef         = N_COEF;
  out->n_coes         = N_COES;
  out->n_disc         = N_DISC;
  out->epsilon        = EPSILON;
  out->earth_radius   = EARTH_RADIUS;
}

int readModelHeader (int *nl, char model[MAX_STRING_LEN],
                     int nd, int Id[nd])
{
  /* Reads velocity model header. */
  /* Open model file */
  FILE *file = fopen (model, "r");

  if (file == NULL) return 1;

  /* Skip first line */
  if (fscanf (file, "%*[^\n]\n") != 0) return 2;
  /* Read number of lines */
  if (fscanf (file, "#N_LINES %d\n", nl) != 1) return 2;
  /* Read discontinuity IDs */
  if (fscanf (file, "#DISC_IDs %d %d %d %d %d %d %d",
              &Id[0], &Id[1], &Id[2], &Id[3],
              &Id[4], &Id[5], &Id[6]) != nd) return 2;

  fclose (file);

  return 0;
}

int readModel (int *ni, int nl,
               char model[MAX_STRING_LEN],
               double rn[nl], double rhon[nl],
               double vpn[nl], double vsn[nl])
{
  /* Reads velocity model. */
  /* Open model file */
  FILE *file = fopen (model, "r");

  if (file == NULL) return 1;

  /* Force C locale so that '.' is used as decimal separator */
  setlocale(LC_NUMERIC, "C");

  /* Skip header lines */
  if (fscanf (file, "%*[^\n]\n") != 0) return 2;
  if (fscanf (file, "%*[^\n]\n") != 0) return 2;
  if (fscanf (file, "%*[^\n]\n") != 0) return 2;
  if (fscanf (file, "%*[^\n]\n") != 0) return 2;

  double depth, rho, vp, vs;

  for (int l = 0; l < nl; l++)
  {
    /* Read line data */
    if (fscanf (file, "%lf %lf %lf %lf",
                &depth, &rho, &vp, &vs) != 4) return 2;

    /* Compute radius */
    rn[l]   = EARTH_RADIUS - depth;
    rhon[l] = rho;
    vpn[l]  = vp;
    vsn[l]  = vs;

    if (l > 0 &&
        fabs (rn[l] - rn[l - 1]) > EPSILON) *ni += 1;
  }

  fclose (file);

  return 0;
}

int checkIO (int rvalue)
{
  /* Checks IO of the model file. */
  switch (rvalue)
  {
    case 1:
      fprintf (stderr, "Error: could not open file!\n");
    break;

    case 2:
      fprintf (stderr, "Error: could not read or write to file!\n");
    break;
  }

  return rvalue;
}