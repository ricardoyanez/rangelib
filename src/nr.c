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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

double nr_polint(double *xa, double *ya, int n, double x, double *dy) {
  double y, c[n], d[n];
  int ns = 1;
  double den, dif, dift, ho, hp, w;
  dif = fabs(x - *xa);
  for ( int i = 0 ; i < n ; i++ ) {
    dift = fabs(x - *(xa+i));
    if ( dift < dif ) {
      ns = i+1;
      dif = dift;
    }
    c[i] = d[i] = *(ya+i);
  }
  y = *(ya+ns-1);
  ns--;
  for ( int m = 0 ; m < n-1 ; m++ ) {
    for ( int i = 0 ; i < n-m-1 ; i++ ) {
      ho = *(xa+i) - x;
      hp = *(xa+i+m+1) - x;
      w = c[i+1] - d[i];
      den = ho - hp;
      if ( den == 0.0 ) {
	fprintf(stderr,"polint: den=0\n");
	exit(EXIT_FAILURE);
      }
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    if ( 2*ns < n-m-1 )
      *dy = c[ns];
    else {
      *dy = d[ns-1];
      ns--;
    }
    y += *dy;
  }
  *dy = fabs(*dy);
  return y;
}

unsigned int nr_locate(double *y, int n, double x) {
  unsigned int jl, jm, ju;
  jl = 0;
  ju = n;
  while ( (ju-jl) > 1 ) {
    jm = (ju + jl) / 2;
    if( (y[n-1]>y[1]) && (x>y[jm]) )
      jl = jm;
    else
      ju = jm;
  }
  return jl;
}

void nr_spline(double *x, double *y, int n, double yp1, double ypn, double *y2) {
  double p, qn, sig, un;
  int m = n-1;
  const int nn = 200;
  double *u = malloc(nn*sizeof(double));
  if ( yp1 > 0.99e30 ) {
    *y2 = 0.0;
    u[0] = 0.0;
  }
  else {
    *y2 = -0.5;
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }
  for ( int k = 1 ; k < n-1 ; k++ ) {
    sig = (x[k] - x[k-1]) / (x[k+1] - x[k-1]);
    p = *(y2+k-1) * sig + 2.0;
    *(y2+k) = (sig - 1.0) / p;
    u[k] = (6.0 * ((y[k+1] - y[k]) / (x[k+1] - x[k]) - 
		   (y[k] - y[k-1]) / (x[k] - x[k-1])) /
	    (x[k+1] - x[k-1]) - sig * u[k-1]) / p;
  }
  if ( ypn > 0.99e30 ) {
    qn = 0.0;
    un = 0.0;
  }
  else {
    qn = 0.5;
    un = (3.0 / (x[m] - x[m-1])) * 
      (ypn - (y[m] - y[m-1]) / (x[m] - x[m-1]));
  }
  *(y2+m) = (un - qn * u[m-1]) / (*(y2+m-1) * qn + 1.0);
  for ( int k = m-1 ; k >= 0 ; k-- ) {
    *(y2+k) *= *(y2+k+1);
    *(y2+k) += u[k];
  }
  free(u);
}

void nr_splint(double *xa, double *ya, double *y2a, int n, double x, double *y) {
  int k, klo, khi;
  double a, b, h;
  klo = 1;
  khi = n;
  while ( khi - klo > 1 ) {
    k = (khi + klo) >> 1;
    if ( xa[k-1] > x ) {
      khi = k;
    }
    else {
      klo = k;
    }
  }
  klo--;
  khi--;
  h = xa[khi] - xa[klo];
  if ( h == 0.0 ) {
    fprintf(stderr,"splint: Bad xa input in splint.\n");
    exit(EXIT_FAILURE);
  }
  a = (xa[khi] - x) / h;
  b = (x - xa[klo]) / h;
  *y = a * ya[klo] + b * ya[khi] + 
    ((pow(a,3) - a) * y2a[klo] + (pow(b,3) - b) * y2a[khi]) * (h*h) / 6.0;
}

void nr_splie2(double *x1a, double *x2a, double **ya,
	       int m, int n, double **y2a) {
  const int nn = 200;
  double *ytmp = malloc(nn*sizeof(double));
  double *y2tmp = malloc(nn*sizeof(double));
  for ( int j = 0 ; j < m ; j++ ) {
    for ( int i = 0 ; i < n ; i++ ) {
      ytmp[i] = *(ya[i]+j);
    }
    nr_spline(x2a,ytmp,n,1.0e30,1.0e30,&y2tmp[0]);
    for ( int i = 0 ; i < n ; i++ ) {
      *(y2a[i]+j) = y2tmp[i];
    }
  }
  free(ytmp);
  free(y2tmp);
}

void nr_splin2(double *x1a, double *x2a, double **ya, double **y2a,
	       int m, int n, double x1, double x2, double *y) {
  const int nn = 200;
  double *ytmp = malloc(nn*sizeof(double));
  double *y2tmp = malloc(nn*sizeof(double));
  double *yytmp = malloc(nn*sizeof(double));
  for ( int j = 0 ; j < m ; j++ ) {
    for ( int i = 0 ; i < n ; i++ ) {
      ytmp[i] = *(ya[i]+j);
      y2tmp[i] = *(y2a[i]+j);
    }
    nr_splint(x2a,ytmp,y2tmp,n,x2,&yytmp[j]);
  }
  nr_spline(x1a,yytmp,m,1.0e30,1.0e30,&y2tmp[0]);
  nr_splint(x1a,yytmp,y2tmp,m,x1,y);
  free(ytmp);
  free(y2tmp);
  free(yytmp);
}

#ifdef __cplusplus
}
#endif
