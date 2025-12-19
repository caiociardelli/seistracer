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

/* Do not change anything in this file unless you know really well what you are doing! */
#ifndef CONSTANTS_H
#define CONSTANTS_H
#define MAX_STRING_LEN_MACRO 200
#define N_BRANCHES_MACRO 4

static const int MAX_STRING_LEN = 200;                            /* Buffer size for short strings */
static const int MAX_PATH_LEN   = 300;                            /* Buffer size for file */
static const int MAX_STREAM_LEN = 100000000;                      /* Buffer size for streams */
static const int MAX_HEAP_DEPTH = 20;                             /* Maximum number of keys in the phase description */
static const int MAX_PHASES     = 5000;                           /* Maximum number of phases per source-receiver pair */
static const int MAX_PATH_NODES = 20000;                          /* Maximum number of nodes in the ray path */
static const int MAX_PAIRS      = 100000;                         /* Maximum number of pairs */
static const int MAX_CODES      = 100000;                         /* Maximum number of permutations */
static const int MAX_ITERATIONS = 10;                             /* Maximum number of iterations to find a solution */
static const int N_TOAG         = 901;                            /* Number of take-off angles */
static const int N_COEF         = 5;                              /* Number of reflection coefficients on the free surface */
static const int N_COES         = 20;                             /* Number of reflection and transmission coefficients at seismic discontinuities */
static const int N_DISC         = 5;                              /* Number of seismic discontinuities */
static const int N_BRANCHES     = 4;                              /* Number of keys "Pu", "Pd", "Su", and "Sd" */
static const int NK             = 4;                              /* The key maps are NK x NK in size */

static const double PI            = 3.14159265358979323846;       /* Pi approximation */
static const double TO_DEGREE     = 180 / 3.14159265358979323846; /* Constant to convert from radians to degrees */
static const double TO_RADIANS    = 3.14159265358979323846 / 180; /* Constant to convert from degrees to radians */
static const double EPSILON       = 1E-15;                        /* Tiny value for doing numerical stuff */
static const double MIN_RADIUS    = 1E-9;                         /* Minimum distance between ray and the Earth's center */
static const double MIN_SLOPE     = 1E-6;                         /* Minimum velocity gradient inside layer */
static const double MIN_PSPH      = 1E-3;                         /* Minimum spherical ray parameter */
static const double MIN_ERROR     = 1E-3;                         /* Minimum error allowed for the source-receiver distance */
static const double MAX_ERROR     = 1E-2;                         /* Maximum error allowed for the source-receiver distance */
static const double EARTH_RADIUS  = 6371.0;                       /* Earth's radius */
static const double MIN_AMPLITUDE = 1E-2;                         /* Minimum allowed amplitude */
static const double MAX_DEPTH     = 700;                          /* Maximum allowed source depth */
#endif
