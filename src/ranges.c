/*
  Author: Ricardo Yanez
  Copyright (c) 2005-2023 Ricardo Yanez <ricardo.yanez@calel.org>

  Calculates energy loss of an ion in an absorber by constructing
  range tables based on either the,

  Northcliffe-Schilling correlations  (E/A < 12 MeV/A)
  (L.C. Northcliffe, R.F. Schilling, Nucl. Data Tables A7, 233, 1970).

  or,

  Hubert-Bimbot-Gauvin correlations (2.5 < E/A < 100 MeV/A)
  (Atomic Data and Nuclear Data Tables, 1990, 46, pp. 1-213).

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

#ifndef _NMAX
#define _NMAX
# define NMAX 4000
#endif

#include "nr.h"

#ifdef __cplusplus
extern "C" {
#endif

void rangetab(int, int, int, int, int, int, double*, double*, int*);

/*
  Calculate energy of ion after passage through an absorber foil.
*/
double passage(int icorr, int zp, int ap, int iabso, int zt, int at,
	       double ein, double t, double *err) {

  double eut, elin, elut, rin, rut, lerr;
  int n;

  int jj, jjj;

  // check correlation
  if ( icorr == 0 && ein/ap > 12.0 ) icorr = 1;  // switch to H-B-G
  if ( icorr == 1 && ein/ap <= 2.5 ) icorr = 0;  // switch to N-S

  double *em = malloc(NMAX*sizeof(double));
  double *r = malloc(NMAX*sizeof(double));

  rangetab(icorr,zp,ap,iabso,zt,at,em,r,&n);
  jjj = n-3;

#ifdef _DEBUG
  FILE *fd;
  if ( icorr == 0 ) {
    fd = fopen("rangetab_ns.dat","w");
  }
  else {
    fd = fopen("rangetab_hbg.dat","w");
  }
  for ( int i = 0 ; i < n ; i++ ) {
    fprintf(fd,"%f\t%f\n",pow(10.0,em[i]),r[i]);
  }
  fclose(fd);
#endif

  elin = log10(ein/ap);
  jj = nr_locate(em,n,elin);
  if ( jj > jjj ) jj = jjj;
  rin = nr_polint(&em[jj],&r[jj],3,elin,&lerr);
  rut = rin - t;
  if ( rut <= 0.0 ) {
    *err = 0.0;
    eut = 0.0;
  }
  else {
    jj = nr_locate(r,n,rut);
    if ( jj > jjj ) jj = jjj;
    elut = nr_polint(&r[jj],&em[jj],3,rut,&lerr);
    *err = fabs(pow(10.0,elut-lerr*3)-pow(10.0,elut+lerr*3))/pow(10.0,elut);
    eut = pow(10.0,elut)*ap;
  }

  // free allocated memory
  free(em);
  free(r);

  return eut;
}

/*
  Calculate incoming energy of ion before passage 
  through an absorber of thickness t.
*/
double egassap(int icorr, int zp, int ap, int iabso, int zt, int at,
	       double t, double eut, double *err) {

  double elut, elin, eaut, rut, rin, lerr;
  int n;

  int jj, jjj;

  double *em = malloc(NMAX*sizeof(double));
  double *r = malloc(NMAX*sizeof(double));

  if ( eut/ap != 0.0 ) {
    if ( icorr == 0 && eut/ap > 12.0 ) icorr = 1;  // switch to H-B-G
  }

  rangetab(icorr,zp,ap,iabso,zt,at,em,r,&n);
  jjj = n-3;

  if ( eut/ap != 0.0 ) {
    elut = log10(eut/ap);
    jj = nr_locate(em,n,elut);
    if (jj > jjj) jj = jjj;
    rut = nr_polint(&em[jj],&r[jj],3,elut,&lerr);
  }
  else {
    rut = 0.0;
  }

  rin = rut + t;
  jj = nr_locate(r,n,rin);
  if ( jj > jjj ) jj = jjj;
  elin = nr_polint(&r[jj],&em[jj],3,rin,&lerr);
  *err = fabs(pow(10.0,elin-lerr*3)-pow(10.0,elin+lerr*3))/pow(10.0,elin);
  eaut = pow(10.0,elin);

  if ( icorr == 0 && eaut > 12.0 ) {
    printf("warning: Hubert-Bimbot-Gauvin correlations should be used in this case.\n");
  }
  if ( icorr == 1 && eaut <= 2.5 ) {
    printf("Warning: Northcliffe-Schilling correlations should be used in this case.\n");
  }

  // free allocated memory
  free(em);
  free(r);

  return eaut*ap;
}

/*
  Calculate absorber thickness for a given energy decrement
*/
double thickn(int icorr, int zp, int ap, int iabso, int zt, int at,
	      double ein, double delen) {

  double elin, elut, rin, rut, rerr;
  int n;

  int jj, jjj;

  double *em = malloc(NMAX*sizeof(double));
  double *r = malloc(NMAX*sizeof(double));

  if ( icorr == 0 && ein/ap > 12.0 ) icorr = 1;  // switch to H-B-G
  if ( icorr == 1 && ein/ap <= 2.5 ) icorr = 0;  // switch to N-S

  rangetab(icorr,zp,ap,iabso,zt,at,em,r,&n);
  jjj = n-3;

  elin = log10(ein/ap);
  jj = nr_locate(em,n,elin);
  if ( jj > jjj ) jj = jjj;
  rin = nr_polint(&em[jj],&r[jj],3,elin,&rerr);
  if ( ein-delen <= 0.0 ) {
    rut = 0.0;
  }
  else {
    elut = log10((ein-delen)/ap);
    jj = nr_locate(em,n,elut);
    if ( jj > jjj ) jj = jjj;
    rut = nr_polint(&em[jj],&r[jj],3,elut,&rerr);
  }

  // free allocated memory
  free(em);
  free(r);

  return rin-rut;
}

/*
  Calculate the range of a projectile
*/
double rangen(int icorr, int zp, int ap, int iabso, int zt, int at,
	      double ein) {

  double rut, elin, rerr;
  int n;

  int jj, jjj;

  double *em = malloc(NMAX*sizeof(double));
  double *r = malloc(NMAX*sizeof(double));

  if ( icorr == 0 && ein/ap > 12.0 ) icorr = 1;  // switch to H-B-G
  if ( icorr == 1 && ein/ap <= 2.5 ) icorr = 0;  // switch to N-S

  rangetab(icorr,zp,ap,iabso,zt,at,em,r,&n);
  jjj = n-3;

  elin = log10(ein/ap);
  jj = nr_locate(em,n,elin);
  if ( jj > jjj ) jj = jjj;
  rut = nr_polint(&em[jj],&r[jj],3,elin,&rerr);

  // free allocated memory
  free(em);
  free(r);

  return rut;
}

#ifdef __cplusplus
}
#endif
