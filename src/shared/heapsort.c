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

static inline void swap (int i, int j, int n, int Ii[n], double T[n])
{
  /* Swaps elements in arrays. */
  int ti = Ii[i]; Ii[i] = Ii[j]; Ii[j] = ti;
  double t = T[i]; T[i] = T[j]; T[j] = t;
}

static void insertIntoHeap (int m, int n, int Ii[n], double T[n])
{
  /* Insert a new element into the max-heap. */
  int p, f = m;

  while (f > 0 && T[p = f / 2] < T[f])
  {
    swap (p, f, n, Ii, T);

    f = p;
  }
}

static void shakeHeap (int m, int n, int Ii[n], double T[n])
{
  /* Shakes array to turn it into a max-heap again after
    having removed the maximum value. */
  int f = 1;

  while (f <= m)
  {
    if (f < m && T[f] < T[f + 1]) f++;
    if (T[f / 2] >= T[f]) break;

    swap (f / 2, f, n, Ii, T);

    f *= 2;
  }
}

void heapSort (int n, int Ii[n], double T[n])
{
  /* Sorts array using the heap sort algorithm. */
  for (int m = 0; m < n; m++)

    insertIntoHeap (m, n, Ii, T);

  for (int m = n - 1; m > 0; m--)
  {
    swap (0, m, n, Ii, T);
    shakeHeap (m - 1, n, Ii, T);
  }
}