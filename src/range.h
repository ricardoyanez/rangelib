/*
 * Copyright (c) 2004-2023 by Ricardo Yanez <ricardo.yanez@calel.org>
 *
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 *
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _RANGE
#define _RANGE
# define NELMAX 10
int nelem;
struct elem {
  int z;
  int a;
  double w;
} absorb[NELMAX];

double passage(int icorr, int zp, int ap, int iabso, int zt, int at,
	       double ein, double t, double *err);

double egassap(int icorr, int zp, int ap, int iabso, int zt, int at,
	       double t, double eut, double *err);

double thickn(int icorr, int zp, int ap, int iabso, int zt, int at,
	      double ein, double delen);

double rangen(int icorr, int zp, int ap, int iabso, int zt, int at,
	      double ein);
#endif

#ifdef __cplusplus
}
#endif
