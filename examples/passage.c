/*
 * Copyright (c) 2004-2023 by Ricardo Yanez <ricardo.yanez@calel.org>
 *
 * Example of how to define compounds and use the rangelib C functions
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

#include <stdio.h>
#include <range.h>

/* Calculate energy loss of alpha particles in air */

int main () {

  int i;
  double dens, t;
  double d;
  double ein, eout, eut, err;

  /* define absorber (air) */
  nelem = 3;
  absorb[0].z = 8;  absorb[0].a = 16; absorb[0].w = 2*16*23.2; /* O2 */
  absorb[1].z = 7;  absorb[1].a = 14; absorb[1].w = 2*14*75.5; /* N2 */
  absorb[2].z = 18; absorb[2].a = 40; absorb[2].w = 1*40*1.3;  /* Ar */

  /* density of air at 1 atm (see README file) */
  dens = 760.0 * 5.38083e-5 * (2*16*23.2 + 2*14*75.5 + 1*40*1.3);

  /* If you want to be more accurate you can use molecular weights instead
     of mass numbers when defining the weights */

  printf("\nCalculating the energy of an alpha particle in air\n\n");

  printf("Initial \t mm of air \t final\n");
  printf("energy  \t           \t energy\n");
  printf("(MeV)   \t           \t [1]    \t [2]\n");
  printf("------- \t --------- \t ------\n");
  ein = 25.0;
  for (i = 0 ; i < 10 ; i++) {
    d = (i+1) / 20.0;
    t = dens * d;
    eout = passage(0,2,4,-1,0,0,ein,t,&err);
    eut = passage(1,2,4,-1,0,0,ein,t,&err);
    printf("%6.1lf \t\t %6.3lf \t %6.2lf \t %6.2lf\n",ein,d*10.0,eout,eut);
  }

  printf("\n[1] Northcliffe-Schilling correlations\n");
  printf("[2] Hubert-Bimbot-Gauvin correlations\n\n");

  return 0;

}
