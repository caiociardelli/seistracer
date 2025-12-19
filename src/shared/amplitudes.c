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
#include <complex.h>
#include "constants.h"
#include "enums.h"
#include "exmath.h"

double complex freeSurface (enum keysF key,
                            double complex a, double complex b,
                            double complex c, double complex d)
{
  /* Computes the reflection amplitude coefficients for a free-surface. */
  double complex asq = csquare (a);

  switch (key)
  {
    case SHfSH:
      return 1;
    break;

    case PfP:
      return (-asq + b) / (asq + b);
    break;

    case PfSV:
      return c * a / (asq + b);
    break;

    case SVfP:
      return d * a / (asq + b);
    break;

    case SVfSV:
      return (asq - b) / (asq + b);
    break;

    default:
      fprintf (stderr, "Unknown coefficient.\n");
      return 0;
    break;
  }
}

double complex solid2Solid (enum keysS key,
                            double p, double psq,
                            double a, double b, double c, double d,
                            double alpha1, double alpha2,
                            double beta1, double beta2,
                            double rho1, double rho2,
                            double complex cosi1, double complex cosi2,
                            double complex cosj1, double complex cosj2,
                            double complex E, double complex F, double complex G,
                            double complex H, double complex D, double complex delta)
{
  /* Computes the reflection and transmission amplitude coefficients
     for a solid-solid interface. */
  double complex k, k1, k2;

  switch (key)
  {
    case SHdSHd:
      return 2 * rho1 * beta1 * cosj1 / delta;
    break;

    case SHdSHu:
      return (rho1 * beta1 * cosj1 - rho2 * beta2 * cosj2) / delta;
    break;

    case SHuSHd:
      return (rho2 * beta2 * cosj2 - rho1 * beta1 * cosj1) / delta;
    break;

    case SHuSHu:
      return 2 * rho2 * beta2 * cosj2 / delta;
    break;

    case PdPd:
      return 2 * rho1 * cosi1 * F / (alpha2 * D);
    break;

    case PdPu:
      k1 = (b * cosi1 / alpha1 - c * cosi2 / alpha2) * F;
      k2 = (a + d * (cosi1 / alpha1) * (cosj2 / beta2)) * H * psq;

      return (k1 - k2) / D;
    break;

    case PuPd:
      k1 = (b * cosi1 / alpha1 - c * cosi2 / alpha2) * F;
      k2 = (a + d * (cosi2 / alpha2) * (cosj1 / beta1)) * G * psq;

      return -1 * (k1 + k2) / D;
    break;

    case PuPu:
      return 2 * rho2 * cosi2 * F / (alpha1 * D);
    break;

    case PdSVd:
      return 2 * rho1 * cosi1 * H * p / (beta2 * D);
    break;

    case PdSVu:
      k = c * d * (cosi2 / alpha2) * (cosj2 / beta2);

      return -2 * cosi1 * p * (a * b + k) / (beta1 * D);
    break;

    case PuSVd:
      k = b * d * (cosi1 / alpha1) * (cosj1 / beta1);

      return 2 * cosi2 * p * (a * c + k) / (beta2 * D);
    break;

    case PuSVu:
      return -2 * rho2 * cosi2 * G * p / (beta1 * D);
    break;

    case SVdSVd:
      return 2 * rho1 * cosj1 * E / (beta2 * D);
    break;

    case SVdSVu:
      k1 = (b * cosj1 / beta1 - c * cosj2 / beta2) * E;
      k2 = (a + d * (cosi2 / alpha2) * (cosj1 / beta1)) * G * psq;

      return -1 * (k1 - k2) / D;
    break;

    case SVuSVd:
      k1 = (b * cosj1 / beta1 - c * cosj2 / beta2) * E;
      k2 = (a + d * (cosi1 / alpha1) * (cosj2 / beta2)) * H * psq;

      return (k1 + k2) / D;
    break;

    case SVuSVu:
      return 2 * rho2 * cosj2 * E / (beta1 * D);
    break;

    case SVdPd:
      return -2 * rho1 * cosj1 * G * p / (alpha2 * D);
    break;

    case SVdPu:
      k = c * d * (cosi2 / alpha2) * (cosj2 / beta2);

      return -2 *cosj1 * p * (a * b + k) / (alpha1 * D);
    break;

    case SVuPd:
      k = b * d * (cosi1 / alpha1) * (cosj1 / beta1);

      return 2 * cosj2 * p * (a * c + k) / (alpha2 * D);
    break;

    case SVuPu:
      return 2 * rho2 * cosj2 * H * p / (alpha1 * D);
    break;

    default:
      fprintf (stderr, "Unknown coefficient.\n");
      return 0;
    break;
  }
}

double complex solid2Liquid (enum keysS key,
                             double p, double a, double c,
                             double alpha1, double alpha2,
                             double beta1, double epsilon,
                             double complex cosi1, double complex cosi2,
                             double complex cosj1,
                             double complex D, double complex E, double complex F,
                             double complex delta)
{
  /* Computes the reflection and transmission amplitude coefficients
     for a solid-liquid interface. */
  switch (key)
  {
    case SHdSHd:
      return 0;
    break;

    case SHdSHu:
      return 1;
    break;

    case PdPd:
      return 2 * cosi1 * a / alpha2 / delta;
    break;

    case PdPu:
      return (D + (E - F) * (cosi2 / alpha2)) / delta;
    break;

    case PuPd:
      return (-D + (E + F) * (cosi2 / alpha2)) / delta;
    break;

    case PuPu:
      return (2 * epsilon / alpha1) * a * cosi2 / delta;
    break;

    case PdSVd:
      return 0;
    break;

    case PdSVu:
      return (-c * cosi1 * a * (cosi2 / alpha2)) / delta;
    break;

    case PuSVd:
      return 0;
    break;

    case PuSVu:
      return -4 * epsilon * beta1 * p * cosi2 * (cosi1 / alpha1) / delta;
    break;

    case SVdSVd:
      return 0;
    break;

    case SVdSVu:
      return (-D + (E - F) * (cosi2 / alpha2)) / delta;
    break;

    case SVdPd:
      return (-c * (beta1 / alpha2) * cosj1 * (cosi1 / alpha1)) / delta;
    break;

    case SVdPu:
      return (beta1 / alpha1) * c * cosj1 * a * (cosi2 / alpha2) / delta;
    break;

    default:
      fprintf (stderr, "Unknown coefficient.\n");
      return 0;
    break;
  }
}

double complex liquid2Solid (enum keysS key,
                             double p, double a, double c,
                             double alpha1, double alpha2,
                             double beta2, double epsilon,
                             double complex cosi1, double complex cosi2,
                             double complex cosj2,
                             double complex D, double complex E, double complex F,
                             double complex delta)
{
  /* Computes the reflection and transmission amplitude coefficients
     for a liquid-solid interface. */
  switch (key)
  {
    case SHuSHu:
      return 0;
    break;

    case SHuSHd:
      return 1;
    break;

    case PdPd:
      return (2 / (alpha2 * epsilon)) * a * cosi1 / delta;
    break;

    case PdPu:
      return (-D + (E + F) * (cosi1 / alpha1)) / delta;
    break;

    case PuPd:
      return (D + (E - F) * (cosi1 / alpha1)) / delta;
    break;

    case PuPu:
      return 2 * cosi2 * a / alpha1 / delta;
    break;

    case PdSVd:
      return -(4 / epsilon) * beta2 * p * cosi1 * (cosi2 / alpha2) / delta;
    break;

    case PdSVu:
      return 0;
    break;

    case PuSVd:
      return (-c * cosi2 * a * (cosi1 / alpha1)) / delta;
    break;

    case PuSVu:
      return 0;
    break;

    case SVuSVd:
      return (-D + (E - F) * (cosi1 / alpha1)) / delta;
    break;

    case SVuSVu:
      return 0;
    break;

    case SVuPd:
      return (beta2 / alpha2) * c * cosj2 * a * (cosi1 / alpha1) / delta;
    break;

    case SVuPu:
      return (-c * (beta2 / alpha1) * cosj2 * (cosi2 / alpha2)) / delta;
    break;

    default:
      fprintf (stderr, "Unknown coefficient.\n");
      return 0;
    break;
  }
}