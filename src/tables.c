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

-----------------------------------------------------------------------------------------------

 TABLES

 USAGE
   ./bin/tables DEPTH MODEL

 EXAMPLE
   ./bin/tables 0 iasp91

 COMMAND-LINE ARGUMENTS
   DEPTH                  - source depth
   MODEL                  - reference 1D seismic velocity model name

 DESCRIPTION
   Creates the tables of ray parameters and distances for a given depth.

----------------------------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "io.h"
#include "permutations.h"

static void helpMenu (void)
{
  /* Prints help menu. */
  char *help_menu = "\n TABLES"

                    "\n\n USAGE"
                    "\n    ./bin/tables DEPTH MODEL"

                    "\n\n EXAMPLE"
                    "\n    ./bin/tables 0 iasp91"

                    "\n\n COMMAND LINE ARGUMENTS"
                    "\n    DEPTH                  - source depth"
                    "\n    MODEL                  - reference 1D seismic velocity model name"

                    "\n\n DESCRIPTION"
                    "\n    Creates the tables of ray parameters and distances for a given depth.\n\n";

  fprintf (stderr, "%s", help_menu);
}

int main (int argc, char *argv[])
{
  if (argc != 3)
  {
    fprintf (stderr, "Error: wrong number of parameters on the command" \
                     " line...\n");
    helpMenu ();

    exit (EXIT_FAILURE);
  }

  /* Parse source depth from command line */
  double depth = atof (argv[1]);

  if (depth > MAX_DEPTH)
  {
    fprintf (stderr, "Error: source depth cannot be larger than %g km!\n",
             MAX_DEPTH); exit (EXIT_FAILURE);
  }

  /* Copy model name from command line */
  char model[MAX_STRING_LEN];

  strncpy (model, argv[2], MAX_STRING_LEN - 1);

  model[MAX_STRING_LEN - 1] = '\0';

  /* Initialize model variables */
  int ni = 0,nl = 0;
  int nd = N_DISC + 2, Id[nd];

  /* Construct model file path */
  char filepath[MAX_PATH_LEN];

  fprintf (stdout, "Reading model %s...\n", model);

  snprintf (filepath, MAX_PATH_LEN, "models/%s.model", model);
  
  /* Read model header */
  readModelHeader (&nl, filepath, nd, Id);

  /* Allocate memory for model arrays */
  double *rn   = malloc (nl * sizeof (double));
  double *rhon = malloc (nl * sizeof (double));
  double *vpn  = malloc (nl * sizeof (double));
  double *vsn  = malloc (nl * sizeof (double));

  if (rn == NULL || rhon == NULL || vpn == NULL || vsn == NULL)
  {
    fprintf (stderr, "Error: could not allocate memory for" \
                     " velocity model...\n"); exit (EXIT_FAILURE);
  }

  /* Read model data */
  if (checkIO (readModel (&ni, nl, filepath,
               rn, rhon, vpn, vsn))) exit (EXIT_FAILURE);

  fprintf (stdout, "Source depth: %.1lf km\n", depth);
  fprintf (stdout, "Computing permutations for each take-off angle...\n");

  /* Initialize id and allocate stream */
  int id; char *stream = malloc (MAX_STREAM_LEN * sizeof (char));

  if (stream == NULL)
  {
    fprintf (stderr, "Error: could not allocate memory for"
                     " stream...\n"); exit (EXIT_FAILURE);
  }

  /* Initialize coefficient arrays */
  double complex cf_FP[N_TOAG][N_COEF];
  double complex cf_FS[N_TOAG][N_COEF];
  double complex cf_SP[N_TOAG][N_DISC][N_COES];
  double complex cf_SS[N_TOAG][N_DISC][N_COES];

  fprintf (stdout, "Precomputing reflection and transmission coefficients...\n");

  /* Precompute all coefficients */
  precomputeAllCoefficients (ni, depth, nl, rn, rhon, vpn, vsn,
                             nd, Id, cf_FP, cf_FS, cf_SP, cf_SS);

  fprintf (stdout, "Creating permutations...\n");

  /* Generate permutations */
  createPermutations (ni, depth, &id, stream,
                      nl, rn, rhon, vpn, vsn,
                      nd, Id, cf_FP, cf_FS, cf_SP, cf_SS);

  /* Free model arrays */
  free (rn);
  free (rhon);
  free (vpn);
  free (vsn);

  fprintf (stdout, "Sorting permutations by source-receiver distance...\n");

  /* Sort permutations */
  if (checkIO (sortByDelta (id, stream, depth))) exit (EXIT_FAILURE);

  /* Free stream memory */
  free (stream);

  fprintf (stdout, "Done!\n");

  return 0;
}