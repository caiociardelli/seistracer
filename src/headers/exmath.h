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

#ifndef EXMATH_H
#define EXMATH_H
static inline double square (double v)
{
  /* Returns the square of a real number. */
  return v * v;
}

static inline double complex csquare (double complex v)
{
  /* Returns the square of a complex number. */
  return v * v;
}

static inline double rad2Degree (double angle)
{
  /* Converts radians to degrees. */
  return angle * TO_DEGREE;
}

static inline double degree2Rad (double v)
{
  /* Converts degrees to radians. */
  return v * TO_RADIANS;
}

static inline double depth2R (double depth)
{
  /* Converts depth to radius. */
  return EARTH_RADIUS - depth;
}
#endif