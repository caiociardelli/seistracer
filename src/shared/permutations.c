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
#include <string.h>
#include <complex.h>
#include <math.h>
#include "constants.h"
#include "enums.h"
#include "structs.h"
#include "io.h"
#include "exmath.h"
#include "source_and_model.h"
#include "amplitudes.h"
#include "ray.h"

static inline void copy (int N, int V1[N], int V2[N])
{
  /* Copies an integer array V1 to another array V2. */
  for (int n = 0; n < N; n++) V2[n] = V1[n];
}

static void preComputeCoeffs (double rs, double vs,
                              struct Surface sfc,
                              struct Discontinuity dsc[N_DISC],
                              double complex cfF[N_TOAG][N_COEF],
                              double complex cfS[N_TOAG][N_DISC][N_COES])
{
  /* Precomputes the reflection and transmission coefficients
     for all discontinuities and take-off angles. */
  for (int tki = 0; tki < N_TOAG; tki++)
  {
    /* Compute take-off angle */
    double take_off = 0.1 * tki;

    /* Compute sin theta */
    double sin_theta = sin (degree2Rad (take_off));

    /* Compute psph */
    double psph = rs * sin_theta / vs;
    double p = psph / sfc.r, psq = square (p);

    /* Get alpha and beta */
    double alpha = sfc.alpha, beta = sfc.beta;

    /* Compute ci and cj */
    double ci = 1 - square (alpha) * psq;
    double cj = 1 - square (beta) * psq;

    /* Compute complex cosines */
    double complex cosi = csqrt (ci);
    double complex cosj = csqrt (cj);

    /* Determine ki and kj */
    double ki = (ci < 0) ? 0 : 1;
    double kj = (cj < 0) ? 0 : 1;

    /* Compute coefficients a, b, c, d */
    double complex a = 1 / square (beta) - 2 * psq;
    double complex b = 4 * psq * cosi * cosj / (alpha * beta);
    double complex c = 4 * (alpha / beta) * p * cosi / alpha;
    double complex d = 4 * (beta / alpha) * p * cosj / beta;

    for (int key = SHfSH; key <= SVfSV; key++)
    {
      /* Compute free surface coefficient */
      double complex cf = freeSurface (key, a, b, c, d);

      if (key == PfSV) cf *= kj; else if (key == SVfP) cf *= ki;

      cfF[tki][key] = cf;
    }

    for (int l = 0; l < N_DISC - 2; l++)
    {
      /* Update p and psq */
      double p = psph / dsc[l].r, psq = square (p);

      /* Get layer properties */
      double alpha1 = dsc[l].alpha1, beta1 = dsc[l].beta1, rho1 = dsc[l].rho1;
      double alpha2 = dsc[l].alpha2, beta2 = dsc[l].beta2, rho2 = dsc[l].rho2;

      /* Compute k1 and k2 */
      double k1 = rho1 * (1 - 2 * square (beta1) * psq);
      double k2 = rho2 * (1 - 2 * square (beta2) * psq);

      /* Compute a, b, c, d */
      double a = k2 - k1;
      double b = k2 + 2 * rho1 * square (beta1) * psq;
      double c = k1 + 2 * rho2 * square (beta2) * psq;
      double d = 2 * (rho2 * square (beta2) - rho1 * square (beta1));

      /* Compute ci and cj for layers */
      double ci1 = 1 - square (alpha1) * psq;
      double ci2 = 1 - square (alpha2) * psq;
      double cj1 = 1 - square (beta1) * psq;
      double cj2 = 1 - square (beta2) * psq;

      /* Compute complex cosines */
      double complex cosi1 = csqrt (ci1);
      double complex cosi2 = csqrt (ci2);
      double complex cosj1 = csqrt (cj1);
      double complex cosj2 = csqrt (cj2);

      /* Determine ki and kj for layers */
      double ki1 = (ci1 < 0) ? 0 : 1;
      double ki2 = (ci2 < 0) ? 0 : 1;
      double kj1 = (cj1 < 0) ? 0 : 1;
      double kj2 = (cj2 < 0) ? 0 : 1;

      /* Compute E, F, G, H */
      double complex E = b * cosi1 / alpha1 + c * cosi2 / alpha2;
      double complex F = b * cosj1 / beta1 + c * cosj2 / beta2;
      double complex G = a - d * (cosi1 / alpha1) * (cosj2 / beta2);
      double complex H = a - d * (cosi2 / alpha2) * (cosj1 / beta1);

      /* Compute D and delta */
      double complex D = E * F + G * H * psq;
      double complex delta = rho1 * beta1 * cosj1 + rho2 * beta2 * cosj2;

      for (int key = SHdSHd; key <= SVuPu; key++)
      {
        /* Compute solid to solid coefficient */
        double complex cf = solid2Solid (key, p, psq, a, b, c, d,
                                         alpha1, alpha2, beta1, beta2, rho1, rho2,
                                         cosi1, cosi2, cosj1, cosj2, E, F, G, H, D, delta);

        /* Adjust coefficient */
        if (key == PdPd)        cf *= ki2;
        else if (key == PdSVu)  cf *= kj1;
        else if (key == PdSVd)  cf *= kj2;
        else if (key == PuPu)   cf *= ki1;
        else if (key == PuSVd)  cf *= kj2;
        else if (key == PuSVu)  cf *= kj1;
        else if (key == SVdPu)  cf *= ki1;
        else if (key == SVdPd)  cf *= ki2;
        else if (key == SVdSVd ||
                 key == SHdSHd) cf *= kj2;
        else if (key == SVuPu)  cf *= ki1;
        else if (key == SVuPd)  cf *= ki2;
        else if (key == SVuSVu ||
                 key == SHuSHu) cf *= kj1;

        cfS[tki][l][key] = cf;
      }
    }

    for (int l = N_DISC - 2; l < N_DISC; l++)
    {
      /* Update p and psq */
      double p = psph / dsc[l].r, psq = square (p);

      /* Get layer properties */
      double alpha1 = dsc[l].alpha1, beta1 = dsc[l].beta1, rho1 = dsc[l].rho1;
      double alpha2 = dsc[l].alpha2, beta2 = dsc[l].beta2, rho2 = dsc[l].rho2;

      /* Compute epsilon */
      double epsilon = rho2 / rho1;

      /* Compute ci and cj for layers */
      double ci1 = 1 - square (alpha1) * psq;
      double ci2 = 1 - square (alpha2) * psq;
      double cj1 = 1 - square (beta1) * psq;
      double cj2 = 1 - square (beta2) * psq;

      /* Compute complex cosines */
      double complex cosi1 = csqrt (ci1);
      double complex cosi2 = csqrt (ci2);
      double complex cosj1 = csqrt (cj1);
      double complex cosj2 = csqrt (cj2);

      /* Determine ki and kj for layers */
      double ki1 = (ci1 < 0) ? 0 : 1;
      double ki2 = (ci2 < 0) ? 0 : 1;
      double kj1 = (cj1 < 0) ? 0 : 1;
      double kj2 = (cj2 < 0) ? 0 : 1;

      if (l == N_DISC - 2)
      {
        /* Compute coefficients for solid to liquid */
        double a = 1 - 2 * square (beta1) * psq;
        double b = 4 * pow (beta1, 4) * psq;
        double c = 4 * beta1 * p;

        double complex D = epsilon * (cosi1 / alpha1);
        double complex E = b * (cosi1 / alpha1) * (cosj1 / beta1);
        double complex F = square (a);

        double complex delta = D + (E + F) * (cosi2 / alpha2);

        for (int key = SHdSHd; key <= SVuPu; key++)
        {
          if (key == SHuSHd ||
              key == SHuSHu ||
              key == SVuSVd ||
              key == SVuSVu ||
              key == SVuPd  ||
              key == SVuPu) continue;

          double complex cf = solid2Liquid (key, p, a, c,
                                            alpha1, alpha2, beta1, epsilon,
                                            cosi1, cosi2, cosj1, D, E, F, delta);

          /* Adjust coefficient */
          if      (key == PdPd)   cf *= ki2;
          else if (key == PdSVu)  cf *= kj1;
          else if (key == PdSVd)  cf *= kj2;
          else if (key == PuPu)   cf *= ki1;
          else if (key == PuSVd)  cf *= kj2;
          else if (key == PuSVu)  cf *= kj1;
          else if (key == SVdPu)  cf *= ki1;
          else if (key == SVdPd)  cf *= ki2;
          else if (key == SVdSVd ||
                   key == SHdSHd) cf *= kj2;

          cfS[tki][l][key] = cf;
        }
      }

      else
      {
        /* Compute coefficients for liquid to solid */
        double a = 1 - 2 * square (beta2) * psq;
        double b = 4 * pow (beta2, 4) * psq;
        double c = 4 * beta2 * p;

        double complex D = (1 / epsilon) * (cosi2 / alpha2);
        double complex E = b * (cosi2 / alpha2) * (cosj2 / beta2);
        double complex F = square (a);

        double complex delta = D + (E + F) * (cosi1 / alpha1);

        for (int key = SHdSHd; key <= SVuPu; key++)
        {
          if (key == SHdSHd ||
              key == SHdSHu ||
              key == SVdSVd ||
              key == SVdSVu ||
              key == SVdPd  ||
              key == SVdPu) continue;

          double complex cf = liquid2Solid (key, p, a, c,
                                            alpha1, alpha2, beta2, epsilon,
                                            cosi1, cosi2, cosj2, D, E, F, delta);

          /* Adjust coefficient */
          if      (key == PdPd)   cf *= ki2;
          else if (key == PdSVu)  cf *= kj1;
          else if (key == PdSVd)  cf *= kj2;
          else if (key == PuPu)   cf *= ki1;
          else if (key == PuSVd)  cf *= kj2;
          else if (key == PuSVu)  cf *= kj1;
          else if (key == SVuPu)  cf *= ki1;
          else if (key == SVuPd)  cf *= ki2;
          else if (key == SVuSVu) cf *= kj1;

          cfS[tki][l][key] = cf;
        }
      }
    }
  }
}

static void phaseName (double psph,
                       int li, int N,
                       int Pm[N],
                       char pnm[MAX_STRING_LEN],
                       double Psph_BP[N_DISC],
                       double Psph_BS[N_DISC])
{
  /* Assigns a name for the phase given a permutation.
     The names mostly follow the IASPEI Standard, but
     new ones had to be invented for phases previously
     unnamed (e.g., a P wave converted to an S when
     transmitted through the Moho, 410-, or 650-km
     discontinuities). */
  /* Copy permutation */
  int Pmc[N];

  for (int i = 0; i < N; i++) Pmc[i] = Pm[i];

  /* Adjust for core phases */
  int pm0 = Pm[0], l = li;

  if (pm0 == Pu || pm0 == Su) l--;

  for (int i = 1; i < N; i++)
  {
    int pm = Pm[i];

    if (l == N_DISC - 1 && pm == Pd)      Pmc[i] = Kd;
    else if (l == N_DISC && pm == Pd)     Pmc[i] = Id;
    else if (l == N_DISC && pm == Sd)     Pmc[i] = Jd;
    else if (l == N_DISC && pm == Pu)     Pmc[i] = Ku;
    else if (l == N_DISC + 1 && pm == Pu) Pmc[i] = Iu;
    else if (l == N_DISC + 1 && pm == Su) Pmc[i] = Ju;

    if (pm == Pu || pm == Su) l--; else l++;
  }

  /* Build phase name */
  char pnmi[MAX_STRING_LEN];

  l = li;

  int i = 0, j = 0, pm = Pmc[i];

  if (pm == Pu || pm == Su) l--;

  while (i < N)
  {
    int dl = (pm == Pu ||
              pm == Su ||
              pm == Ku ||
              pm == Iu ||
              pm == Ju) ? -1 : 1;

    while (i < N && Pmc[i] == pm)
    {
      if (i > 0) {l += dl;} i++;
    }

    /* Append branch to name */
    switch (Pmc[i - 1])
    {
      case Pu:
        pnmi[j] = 'P'; pnmi[j + 1] = 'u'; j += 2;
      break;

      case Pd:
        pnmi[j] = 'P'; pnmi[j + 1] = 'd'; j += 2;
      break;

      case Su:
        pnmi[j] = 'S'; pnmi[j + 1] = 'u'; j += 2;
      break;

      case Sd:
        pnmi[j] = 'S'; pnmi[j + 1] = 'd'; j += 2;
      break;

      case Ku:
        pnmi[j] = 'K'; pnmi[j + 1] = 'u'; j += 2;
      break;

      case Kd:
        pnmi[j] = 'K'; pnmi[j + 1] = 'd'; j += 2;
      break;

      case Iu:
        pnmi[j] = 'I'; pnmi[j + 1] = 'u'; j += 2;
      break;

      case Id:
        pnmi[j] = 'I'; pnmi[j + 1] = 'd'; j += 2;
      break;

      case Ju:
        pnmi[j] = 'J'; pnmi[j + 1] = 'u'; j += 2;
      break;

      case Jd:
        pnmi[j] = 'J'; pnmi[j + 1] = 'd'; j += 2;
      break;
    }

    if (i == N) break;

    /* Append layer number */
    sprintf (&pnmi[j], "%d", l);

    j++; pm = Pmc[i];
  }

  pnmi[j] = '\0';

  /* Finalize phase name */
  int k = 0;

  if (pnmi[1] == 'u')

    pnm[k] = (pnmi[0] == 'P') ? 'p' : 's';

  else

    pnm[k] = (pnmi[0] == 'P') ? 'P' : 'S';

  k++;

  for (int m = 2; m < j - 2; m += 3)
  {
    int l = atoi (&pnmi[m]);

    if (k > 0 &&
        pnmi[m - 1] == 'd' &&
        pnmi[m + 2] == 'u' &&
        pnmi[m - 2] == pnmi[m + 1])
    {
      char pm = pnmi[m + 1];

      if (l == N_DISC + 1 ||
          (((pm == 'P' ||
             pm == 'K') && psph > Psph_BP[l - 1]) ||
            (pm == 'S'  && psph > Psph_BS[l - 1]))) continue;
    }

    /* Append discontinuity labels */
    switch (l)
    {
      case 1:
        pnm[k] = 'm'; k++;
      break;

      case 2:
        pnm[k] = '4'; k++;
        pnm[k] = '1'; k++;
        pnm[k] = '0'; k++;
      break;

      case 3:
        pnm[k] = '6'; k++;
        pnm[k] = '6'; k++;
        pnm[k] = '0'; k++;
      break;

      case 4:
        if (pnmi[m - 1] == 'd' &&
            pnmi[m + 2] == 'u') {pnm[k] = 'c'; k++;}
      break;

      case 5:
        if (pnmi[m - 1] == 'd' &&
            pnmi[m + 2] == 'u') {pnm[k] = 'i'; k++;}
      break;
    }

    if (l == 1 || l == 2 || l == 3)
    {
      if (pnmi[m - 1] == pnmi[m + 2] &&
          pnmi[m - 2] != pnmi[m + 1])
      {
        pnm[k] = 't'; k++;
        pnm[k] = (pnmi[m - 1] == 'u') ? '+' : '-';
      }

      else if (pnmi[m - 1] == 'd' &&
               pnmi[m + 2] == 'u') pnm[k] = '+';

      else pnm[k] = '-';

      k++;
    }

    if (k > 0) {pnm[k] = pnmi[m + 1]; k++;}
  }

  pnm[k] = '\0';
}

static void initializeTree (struct Node **root, int l, int pm)
{
  /* Initializes the tree to store the rays paths. */
  *root = malloc (sizeof (struct Node));

  /* Set initial node properties */
  (*root)->p  =  0;
  (*root)->l  =  l;
  (*root)->k  =  4;
  (*root)->pm = pm;

  (*root)->is_SH = (pm == Su ||
                    pm == Sd) ? true : false;

  (*root)->c_AV = 1.0;
  (*root)->c_AH = 1.0;

  (*root)->up = NULL;

  for (int b = 0; b < N_BRANCHES; b++)

    (*root)->P[b] = NULL;
}

static void buildTree (struct Node *root, int *M, int tki,
                       int Map_FH[NK][NK], int Map_FV[NK][NK],
                       int Map_SH[NK][NK], int Map_SV[NK][NK],
                       double psph,
                       double PsphBP[N_DISC], double PsphBS[N_DISC],
                       double complex cfF[N_TOAG][N_COEF],
                       double complex cfS[N_TOAG][N_DISC][N_COES])
{
  /* Builds the tree. */
  struct Node *node = root;

  do
  {
    do
    {
      node->k--; int l = node->l;

      if (l >= 0 && l <= 6)
      {
        int pm1 = node->pm;
        int pm2 = node->k;

        bool is_Not_Turning_Layer = true;
        bool is_Valid_Permutation = true;

        if (pm1 == Pd || pm1 == Sd)
        {
          double psphB = (pm1 == Pu ||
                          pm1 == Pd) ? PsphBP[l - 1]
                                     : PsphBS[l - 1];

          if (l == N_DISC + 1 || psph > psphB)
          {
            if ((pm1 == Pd && pm2 == Pu) ||
                (pm1 == Sd && pm2 == Su))

              is_Not_Turning_Layer = false;

            else is_Valid_Permutation = false;
          }
        }

        if (is_Valid_Permutation)
        {
          double c_AV = node->c_AV;
          double c_AH = node->c_AH;

          if (is_Not_Turning_Layer)
          {
            if (l == 0)
            {
              /* Get free surface keys */
              int keyv = Map_FV[pm1][pm2];
              int keyh = Map_FH[pm1][pm2];

              c_AV *= cabs (cfF[tki][keyv]);
              c_AH *= cabs (cfF[tki][keyh]);
            }

            else
            {
              /* Get internal keys */
              int keyv = Map_SV[pm1][pm2];
              int keyh = Map_SH[pm1][pm2];

              c_AV *= cabs (cfS[tki][l - 1][keyv]);
              c_AH *= cabs (cfS[tki][l - 1][keyh]);
            }
          }

          /* Determine if SH */
          bool is_SH = (node->is_SH &&
                        (pm2 == Su ||
                         pm2 == Sd)) ? true : false;

          double cA = is_SH ? c_AH : c_AV;

          if (cA >= MIN_AMPLITUDE)
          {
            /* Allocate new node */
            struct Node *new = malloc (sizeof (struct Node));

            /* Update layer */
            if (node->k < 2) l--; else l++;

            new->l  = l;
            new->p  = node->p + 1;
            new->pm = node->k;
            new->k  = (new->p == MAX_HEAP_DEPTH - 1) ? 0 : 4;

            new->is_SH = is_SH;

            new->c_AV = c_AV;
            new->c_AH = c_AH;

            new->up = node;

            for (int k = 0; k < N_BRANCHES; k++)

              new->P[k] = NULL;

            /* Add new node */
            node->P[node->k] = new;
            node = new;

            if (l == 0) *M += 1;
          }
        }
      }

      while (node->k == 0 &&
             node->up != NULL)

        node = node->up;

    } while (node->up != NULL);

  } while (root->k > 0);
}

static void treeToArray (struct Node *root, int *m,
                         int M, int Lp[M],
                         int Pm[M][MAX_HEAP_DEPTH])
{
  /* Converts the tree to an array. */
  struct Node *node = root;

  int Pmt[MAX_HEAP_DEPTH];

  while (root->k < N_BRANCHES)
  {
    Pmt[root->p] = root->pm;

    while (node->k < N_BRANCHES)
    {
      if (node->P[node->k] != NULL &&
          node->P[node->k]->k < N_BRANCHES)
      {
        /* Traverse down */
        node = node->P[node->k];

        Pmt[node->p] = node->pm;
      }

      else node->k++;
    }

    if (node->l == 0)
    {
      /* Store path */
      Lp[*m] = node->p + 1;

      copy (Lp[*m], Pmt, Pm[*m]);

     *m += 1;
    }

    if (node->up != NULL)
    {
      node = node->up;

      /* Free child */
      free (node->P[node->k]);
      node->P[node->k] = NULL;
    }
  }

  /* Free root */
  free (root); root = NULL;
}

void precomputeAllCoefficients (int ni, double depth,
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
{
  /* Precomputes all coefficients. */
  /* Compute source radius */
  double rs = depth2R (depth);

  double vs_P, vs_S;

  /* Get source velocities */
  sourceVelocity (rs, nl, rn, vpn, vsn, &vs_P, &vs_S);

  /* Create model structure */
  struct Model Md[ni];

  arraysToStructure (ni, Md, nl, rn, vpn, vsn);

  double Psph_BP[N_DISC], Psph_BS[N_DISC];

  struct Surface sfc;
  struct Discontinuity dsc[N_DISC];

  /* Set interfaces */
  setInterfaces (nl, rn, rhon, vpn, vsn, nd, Id,
                 Psph_BP, Psph_BS, &sfc, dsc);

  /* Compute for P */
  preComputeCoeffs (rs, vs_P, sfc, dsc, cf_FP, cf_SP);
  /* Compute for S */
  preComputeCoeffs (rs, vs_S, sfc, dsc, cf_FS, cf_SS);
}

void createPermutations (int ni, double depth, int *id,
                         char stream[MAX_STREAM_LEN],
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
{
  /* Creates all possible permutations of ray paths. */
  /* Compute source radius and delta s */
  double rs = depth2R (depth);
  double ds = 0.5 * PI;

  /* Define keys tags */
  char keys_tags[] = {'A', 'B', 'C', 'D'};

  /* Define mapping arrays */
  int Map_FH[NK][NK];

  Map_FH[Pd][Pd] = -1;     Map_FH[Pu][Pd] = -1;
  Map_FH[Pd][Pu] = -1;     Map_FH[Pu][Pu] = -1;
  Map_FH[Pd][Sd] = -1;     Map_FH[Pu][Sd] = -1;
  Map_FH[Pd][Su] = -1;     Map_FH[Pu][Su] = -1;

  Map_FH[Sd][Pd] = -1;     Map_FH[Su][Pd] = -1;
  Map_FH[Sd][Pu] = -1;     Map_FH[Su][Pu] = -1;
  Map_FH[Sd][Sd] = -1;     Map_FH[Su][Sd] = SHfSH;
  Map_FH[Sd][Su] = -1;     Map_FH[Su][Su] = -1;

  int Map_FV[NK][NK];

  Map_FV[Pd][Pd] = -1;     Map_FV[Pu][Pd] = PfP;
  Map_FV[Pd][Pu] = -1;     Map_FV[Pu][Pu] = -1;
  Map_FV[Pd][Sd] = -1;     Map_FV[Pu][Sd] = PfSV;
  Map_FV[Pd][Su] = -1;     Map_FV[Pu][Su] = -1;

  Map_FV[Sd][Pd] = -1;     Map_FV[Su][Pd] = SVfP;
  Map_FV[Sd][Pu] = -1;     Map_FV[Su][Pu] = -1;
  Map_FV[Sd][Sd] = -1;     Map_FV[Su][Sd] = SVfSV;
  Map_FV[Sd][Su] = -1;     Map_FV[Su][Su] = -1;

  int Map_SH[NK][NK];

  Map_SH[Pd][Pd] = -1;     Map_SH[Pu][Pd] = -1;
  Map_SH[Pd][Pu] = -1;     Map_SH[Pu][Pu] = -1;
  Map_SH[Pd][Sd] = -1;     Map_SH[Pu][Sd] = -1;
  Map_SH[Pd][Su] = -1;     Map_SH[Pu][Su] = -1;

  Map_SH[Sd][Pd] = -1;     Map_SH[Su][Pd] = -1;
  Map_SH[Sd][Pu] = -1;     Map_SH[Su][Pu] = -1;
  Map_SH[Sd][Sd] = SHdSHd; Map_SH[Su][Sd] = SHuSHd;
  Map_SH[Sd][Su] = SHdSHu; Map_SH[Su][Su] = SHuSHu;

  int Map_SV[NK][NK];

  Map_SV[Pd][Pd] = PdPd;   Map_SV[Pu][Pd] = PuPd;
  Map_SV[Pd][Pu] = PdPu;   Map_SV[Pu][Pu] = PuPu;
  Map_SV[Pd][Sd] = PdSVd;  Map_SV[Pu][Sd] = PuSVd;
  Map_SV[Pd][Su] = PdSVu;  Map_SV[Pu][Su] = PuSVu;

  Map_SV[Sd][Pd] = SVdPd;  Map_SV[Su][Pd] = SVuPd;
  Map_SV[Sd][Pu] = SVdPu;  Map_SV[Su][Pu] = SVuPu;
  Map_SV[Sd][Sd] = SVdSVd; Map_SV[Su][Sd] = SVuSVd;
  Map_SV[Sd][Su] = SVdSVu; Map_SV[Su][Su] = SVuSVu;

  double vs_P, vs_S;

  /* Get source velocities */
  sourceVelocity (rs, nl, rn, vpn, vsn, &vs_P, &vs_S);

  /* Create model structure */
  struct Model Md[ni];

  arraysToStructure (ni, Md, nl, rn, vpn, vsn);

  double Psph_BP[N_DISC], Psph_BS[N_DISC];

  struct Surface sfc;
  struct Discontinuity dsc[N_DISC];

  /* Set interfaces */
  setInterfaces (nl, rn, rhon, vpn, vsn, nd, Id,
                 Psph_BP, Psph_BS, &sfc, dsc);

  *id = 0;

  for (int tki = 0; tki < N_TOAG; tki++)
  {
    /* Compute take-off angle */
    double take_off = 0.1 * tki;

    if (tki % 100 == 0)

      fprintf (stdout, "Take-off angle: %4.1lf degrees\n", take_off);

    /* Compute sin theta */
    double sin_theta = sin (degree2Rad (take_off));

    /* Compute psph for P and S */
    double psph_P = rs * sin_theta / vs_P;
    double psph_S = rs * sin_theta / vs_S;

    /* Initialize tree roots */
    struct Node *root_Pu = NULL;
    struct Node *root_Su = NULL;
    struct Node *root_Pd = NULL;
    struct Node *root_Sd = NULL;

    int M = 0, li = setLayer (rs, nl, rn, nd, Id), mi = 0;

    while (Md[mi].rb >= rs) mi++;

    if (rs < EARTH_RADIUS)
    {
      /* Initialize upward trees */
      initializeTree (&root_Pu, li - 1, Pu);
      initializeTree (&root_Su, li - 1, Su);

      if (li == 1) M += 2;

      /* Build P upward tree */
      buildTree (root_Pu, &M, tki,
                 Map_FH, Map_FV,
                 Map_SH, Map_SV,
                 psph_P,
                 Psph_BP, Psph_BS,
                 cf_FP, cf_SP);
      /* Build S upward tree */
      buildTree (root_Su, &M, tki,
                 Map_FH, Map_FV,
                 Map_SH, Map_SV,
                 psph_S,
                 Psph_BP, Psph_BS,
                 cf_FS, cf_SS);
      
      if (root_Pu == NULL || root_Su == NULL)
      {
        fprintf (stderr, "Error: could not allocate memory for"
                         " trees...\n"); exit (EXIT_FAILURE);
      }
    }

    /* Initialize downward trees */
    initializeTree (&root_Pd, li, Pd);
    initializeTree (&root_Sd, li, Sd);

    /* Build P downward tree */
    buildTree (root_Pd, &M, tki,
               Map_FH, Map_FV,
               Map_SH, Map_SV,
               psph_P,
               Psph_BP, Psph_BS,
               cf_FP, cf_SP);
    /* Build S downward tree */
    buildTree (root_Sd, &M, tki,
               Map_FH, Map_FV,
               Map_SH, Map_SV,
               psph_S,
               Psph_BP, Psph_BS,
               cf_FS, cf_SS);

    if (root_Pd == NULL || root_Sd == NULL)
    {
      fprintf (stderr, "Error: could not allocate memory for"
                       " trees...\n"); exit (EXIT_FAILURE);
    }

    int m = 0, Lp[M], Pm[M][MAX_HEAP_DEPTH];

    if (rs < EARTH_RADIUS)
    {
      /* Convert P upward tree to array */
      treeToArray (root_Pu, &m, M, Lp, Pm);
      /* Convert S upward tree to array */
      treeToArray (root_Su, &m, M, Lp, Pm);
    }

    /* Convert P downward tree to array */
    treeToArray (root_Pd, &m, M, Lp, Pm);
    /* Convert S downward tree to array */
    treeToArray (root_Sd, &m, M, Lp, Pm);

    char pnm[MAX_STRING_LEN];

    /* Append number of intervals */
    *id += sprintf (&stream[*id], "# %d\n", ni);

    for (int m = 0; m < M; m++)
    {
      /* Select psph */
      double psph = (Pm[m][0] == Pd ||
                     Pm[m][0] == Pu) ? psph_P : psph_S;

      /* Compute delta */
      double delta = computeDelta (li, mi, ni, Md, Lp[m], Pm[m],
                                   nl, rn, nd, Id, psph, rs, ds);

      /* Append permutation data */
      *id += sprintf (&stream[*id], "%d ", Lp[m]);

      for (int n = 0; n < Lp[m]; n++)

        *id += sprintf (&stream[*id], "%c", keys_tags[Pm[m][n]]);

      /* Get phase name */
      phaseName (psph, li, Lp[m], Pm[m], pnm, Psph_BP, Psph_BS);

      *id += sprintf (&stream[*id], " %s %lf %lf\n", pnm, delta, psph);

      if (*id > MAX_STREAM_LEN)
      {
        fprintf (stderr, "Error: stream length is larger than MAX_STREAM_LEN!\n");

        exit (EXIT_FAILURE);
      }
    }
  }
}

static bool isInPhasesList (int c, int *index,
                            struct Phases pairs[MAX_PAIRS],
                            char phase[MAX_STRING_LEN])
{
  /* Checks whether a phase is already on the list. */
  *index = 0;

  while (*index < c && strcmp (phase,
                               pairs[*index].phase)) *index += 1;

  if (c == *index) return false;

  return true;
}

static bool isNotInCodesList (int d, int *index,
                              char codes[MAX_CODES][MAX_HEAP_DEPTH + 1],
                              char code[MAX_STRING_LEN])
{
 /* Checks whether a code is not on the list. */
  *index = 0;

  while (*index < d && strcmp (code, codes[*index])) *index += 1;

  if (d == *index) return true;

  return false;
}

int sortByDelta (int id, char stream[MAX_STREAM_LEN], double depth)
{
  /* Sorts permutations by their corresponding delta (i.e., the great
     circle arc between source and receiver). */
  char line[MAX_STRING_LEN];
  char code[MAX_STRING_LEN];
  char phase[MAX_STRING_LEN];

  /* Define file names */
  char filename1[MAX_STRING_LEN];
  char filename2[MAX_STRING_LEN];
  char filename3[MAX_STRING_LEN];

  struct Phases pairs[MAX_PAIRS];

  int dp = (int) round (depth);

  if (sprintf (filename1, "tables/codes_%d.list", dp)  < 18) return 1;
  if (sprintf (filename2, "tables/phases_%d.list", dp) < 19) return 1;
  if (sprintf (filename3, "tables/table_%d.bin", dp)   < 18) return 1;

  /* Open output files */
  FILE *file1 = fopen (filename1, "w");
  FILE *file2 = fopen (filename2, "w");
  FILE *file3 = fopen (filename3, "wb");

  if (file1 == NULL) return 2;
  if (file2 == NULL) return 2;
  if (file3 == NULL) return 2;

  /* Initialize pairs */
  for (int i = 0; i < MAX_PAIRS; i++) pairs[i].n = 0;

  int N = 0, c = 0, index;

  for (int j = 0; j < id; j++)
  {
    int k = 0;

    while (stream[j] != '\n') line[k++] = stream[j++];

    line[k] = '\0';

    if (line[0] != '#')
    {
      int np;

      double psph, delta;

      /* Parse line data */
      sscanf (line, "%d %s %s %lf %lf", &np, code, phase,
                                        &delta, &psph);

      if (isInPhasesList (c, &index, pairs, phase)) pairs[index].n++;

      else
      {
        /* Add new phase */
        strcpy (pairs[c].phase, phase); pairs[c++].n++;

        if (c == MAX_PAIRS)
        {
          fprintf (stderr, "Error: number of pairs is larger than MAX_PAIRS!\n");

          exit (EXIT_FAILURE);
        }
      }

      N++;
    }
  }

  /* Compute indices */
  int ii[c]; ii[0] = 0;

  for (int i = 1; i < c; i++)

    ii[i] = ii[i - 1] + pairs[i - 1].n;

  int d = 0, max_code_len = MAX_HEAP_DEPTH + 1;

  /* Allocate codes array */
  char (*codes)[max_code_len] = malloc (sizeof (char[MAX_CODES][max_code_len]));

  if (codes == NULL)
  {
    fprintf (stderr, "Error: could not allocate memory for"
                     " codes...\n"); exit (EXIT_FAILURE);
  }

  /* Allocate blocks */
  struct Block *blocks = malloc (sizeof (struct Block[N]));

  if (blocks == NULL)
  {
    fprintf (stderr, "Error: could not allocate memory for"
                     " blocks...\n"); exit (EXIT_FAILURE);
  }

  for (int j = 0; j < id; j++)
  {
    int k = 0;

    while (stream[j] != '\n') line[k++] = stream[j++];

    line[k] = '\0';

    char hashtag; int nl;

    if (line[0] == '#')

      sscanf (line, "%c %d", &hashtag, &nl);

    else
    {
      int np; double delta, psph;

      sscanf (line, "%d %s %s %lf %lf", &np, code, phase,
                                        &delta, &psph);

      if (isNotInCodesList (d, &index, codes, code))

        strcpy (codes[d++], code);

      if (d == MAX_CODES)
      {
        fprintf (stderr, "Error: number of codes is larger than MAX_CODES!\n");

        exit (EXIT_FAILURE);
      }

      int i = 0;

      while (i < c)
      {
        if (!strcmp (phase, pairs[i].phase)) break;

        i++;
      }

      /* Store block data */
      blocks[ii[i]].np    = (unsigned short) np;
      blocks[ii[i]].index = (int) index;
      blocks[ii[i]].delta = (float) delta;
      blocks[ii[i]].psph  = (float) psph;

      ii[i] += 1;
    }
  }

  /* Write number of codes */
  fprintf (file1, "%d\n", d);

  for (int i = 0; i < d; i++)

    fprintf (file1, "%s\n", codes[i]);

  free (codes);

  /* Compute Nv and jj */
  int Nv = 0, jj[N], kn = 1, kt = 1;

  for (int i = 0, j = 0; i < N; i++, j++)
  {
    jj[Nv] = 0;

    if (j == pairs[kt - 1].n)
    {
      if (pairs[kt].n > 1) kn++;

      j = 0; kt++;
    }

    if (pairs[kt - 1].n > 1)

      jj[Nv++] = i;
  }

  /* Write number of phases */
  fprintf (file2, "%d\n", kn);

  for (int k = 0; k < kt; k++)

    if (pairs[k].n > 1) fprintf (file2, "%d %s\n",
                                 pairs[k].n, pairs[k].phase);

  /* Allocate selected blocks */
  struct Block *selected_blocks = malloc (sizeof (struct Block[N]));

  if (selected_blocks == NULL)
  {
    fprintf (stderr, "Error: could not allocate memory for"
                     " selected_blocks...\n"); exit (EXIT_FAILURE);
  }

  for (int j = 0; j < Nv; j++)

    selected_blocks[j] = blocks[jj[j]];

  free (blocks);

  /* Write number of blocks */
  fwrite (&Nv, sizeof (int), 1, file3);
  /* Write blocks */
  fwrite (selected_blocks, sizeof (struct Block), Nv, file3);

  free (selected_blocks);

  fclose (file1);
  fclose (file2);
  fclose (file3);

  return 0;
}