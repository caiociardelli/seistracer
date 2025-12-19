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

#ifndef STRUCTS_H
#define STRUCTS_H

#include <stdbool.h>

struct SeisTracerConstants
{
  /* SeisTracer constants */
  int    max_string_len;
  int    max_phases;
  int    max_path_nodes;
  int    n_toag;
  int    n_coef;
  int    n_coes;
  int    n_disc;
  double epsilon;
  double earth_radius;
};

struct Model
{
  /* Cell to store model parameters. */
  double rb;
  double rt;
  double vpb;
  double vpt;
  double vsb;
  double vst;
  double avp;
  double bvp;
  double avs;
  double bvs;
};

struct Amplitude
{
  /* Cell to store amplitude values. */
  double v;
  double h;
};

struct Surface
{
  /* Cell to store model parameters on
     the surface of the Earth. */
  double r;
  double rho;
  double alpha;
  double beta;
};

struct Discontinuity
{
  /* Cell to store model parameters at
     the internal discontinuities of
     the Earth. */
  double r;
  double rho1;
  double rho2;
  double alpha1;
  double alpha2;
  double beta1;
  double beta2;
  double pb_alpha;
  double pb_beta;
};

struct Phases
{
  /* Cell to store phase codes. */
  char phase[MAX_STRING_LEN_MACRO];

  int n;
};

struct Node
{
  /* Node for the tree. */
  int p;
  int l;
  int k;
  int pm;

  bool is_SH;

  double c_AV;
  double c_AH;

  struct Node *up;
  struct Node *P[N_BRANCHES_MACRO];
};

struct Block
{
  /* Block of information for the tables. */
  unsigned short np;

  int index;

  float delta;
  float psph;
};
#endif