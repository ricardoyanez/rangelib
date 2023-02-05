/*
  Author: Ricardo Yanez

  Copyright (c) 2004-2023 Ricardo Yanez <ricardo.yanez@calel.org>

  Numerical functions for rangelib: Numerical Recipes in Fortran ported to C.

  License:

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

*/

double nr_polint(double *xa, double *ya, int n, double x, double *dy);
unsigned int nr_locate(double *y, int n, double x);
void nr_splie2(double *x1a, double *x2a, double **ya,
	       int m, int n, double **y2a);
void nr_splin2(double *x1a, double *x2a, double **ya, double **y2a,
	       int m, int n, double x1, double x2, double *y);
