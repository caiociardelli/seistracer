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

#ifndef ENUMS_H
#define ENUMS_H
/* Direction keys */
enum keysB {U, D};

/* Wave type keys */
enum keysK {P, S};

/* Permutation keys */
enum keysP {Pu, Su,
            Pd, Sd,
            Ku, Kd,
            Iu, Ju,
            Id, Jd,
            Skip};

/* Free surface coefficient keys */
enum keysF {SHfSH, PfP, PfSV, SVfP, SVfSV};

/* Internal coefficient keys */
enum keysS {SHdSHd, SHdSHu, SHuSHd, SHuSHu,
            PdPd, PdPu, PuPd, PuPu,
            PdSVd, PdSVu, PuSVd, PuSVu,
            SVdSVd, SVdSVu, SVuSVd, SVuSVu,
            SVdPd, SVdPu, SVuPd, SVuPu};
#endif