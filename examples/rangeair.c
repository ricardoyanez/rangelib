/*
 * Copyright (c) 2022-2023 by Ricardo Yanez <ricardo.yanez@calel.org>
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

/* Calculate the range of protons in air */

int main () {

  int i;
  double dens;
  double ein, rout, rut;

  /* define absorber (air) */
  nelem = 3;
  absorb[0].z = 8;  absorb[0].a = 16; absorb[0].w = 2*16*23.2; /* O2 */
  absorb[1].z = 7;  absorb[1].a = 14; absorb[1].w = 2*14*75.5; /* N2 */
  absorb[2].z = 18; absorb[2].a = 40; absorb[2].w = 1*40*1.3;  /* Ar */

  /* density of air at 1 atm (see README file) */
  dens = 760.0 * 5.38083e-5 * (2*16*23.2 + 2*14*75.5 + 1*40*1.3);

  /* If you want to be more accurate you can use molecular weights instead
     of mass numbers when defining the weights */

  printf("\nCalculating the range of protons in air\n\n");

  printf("Initial \t mm of air \t mm of air\n");
  printf("energy  \t [1] \t\t [2]\n");
  printf("(MeV)   \t \n");
  printf("------- \t --------- \t --------- \n");
  ein = 0.0;
  for (i = 0 ; i < 20 ; i++) {
    ein += 1.0;
    rout = rangen(0,1,1,-1,0,0,ein) / dens;
    rut = rangen(1,1,1,-1,0,0,ein) / dens;
    printf("%6.1lf \t\t %6.2lf\t\t %6.2lf\n",ein,rout*10.0,rut*10.0);
  }

  printf("\n[1] Northcliffe-Schilling correlations\n");
  printf("[2] Hubert-Bimbot-Gauvin correlations\n\n");

  return 0;

}
