/*
  Author: Ricardo Yanez

  Copyright (c) 2004-2023 Ricardo Yanez <ricardo.yanez@calel.org>

  Program for -dE/dx and range calculations for ions in various
  absorber foils.

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
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#ifndef _NMAX
#define _NMAX
# define NMAX 4000
#endif

#include "range.h"
#include "nr.h"

// major version number
static const char ver[] = "0.2.0"; 

struct elem cmpnd[NELMAX];
int numel;

void rangetab(int, int, int, int, int, int, double*, double*, int*);
void dedxtab(int, int, int, int, int, int, double, double*, double*);

/*
 * Calculate range for a given initial energy
 */
double rangestop(int icorr, int zp, int ap, int iabso, int zt, int at, double ein) {

  double delr, elin, err;
  int n;

  int jj, jjj;

  double *em = malloc(NMAX*sizeof(double));
  double *r = malloc(NMAX*sizeof(double));

  if ( icorr == 0 && ein/ap > 12.0 ) {
    icorr = 1;  // switch to H-B-G
  }
  if ( icorr == 1 && ein/ap <= 2.5 ) {
    icorr = 0;  // switch to N-S
  }

  rangetab(icorr,zp,ap,iabso,zt,at,&em[0],&r[0],&n);
  jjj = n-3;

  elin = log10(ein/ap);
  jj = nr_locate(em,n,elin);
  if ( jj > jjj ) jj = jjj;
  delr = nr_polint(&em[jj],&r[jj],3,elin,&err);

  // free allocated memory
  free(em);
  free(r);

  return delr;
}

/*
 * Calculate initial energy from a given range
 */
double rangein(int icorr, int zp, int ap, int iabso, int zt, int at, double t) {

  double ein, elin, err;
  int n;

  int jj, jjj;

  double *em = malloc(NMAX*sizeof(double));
  double *r = malloc(NMAX*sizeof(double));

  rangetab(icorr,zp,ap,iabso,zt,at,&em[0],&r[0],&n);
  jjj = n-3;

  jj = nr_locate(r,n,t);
  if ( jj > jjj ) jj = jjj;
  elin = nr_polint(&r[jj],&em[jj],3,t,&err);
  ein = pow(10,elin) * ap;

  // free allocated memory
  free(em);
  free(r);

  return ein;
}

void disp_copy() {
  printf("Copyright(c) 2004-2023 Ricardo Yanez <ricardo.yanez@calel.org>\n");
}

void disp_header() {
  printf("Welcome to the -dE/dx and range program version %s\n",ver);
  disp_copy();
}

void disp_help() {
  printf("\nUsage: range [OPTION]\n");
  printf("Options:\n");
  printf("  -n, --ns         Use Northcliffe-Schilling correlations (default)\n");
  printf("  -h, --hbg        Use Hubert-Bimbot-Gauvin correlations\n");
  printf("  -l, --list       List available compounds and exit\n");
  printf("  -v, --version    Display rangelib version number and exit\n");
  printf("      --help       Display this help and exit\n\n");
}

void list_compounds() {
  printf("\n\t  -1: User-defined compound\n");
  printf("\n\tSolids:\n");
  printf("\t      Compound\t\t\t\tDensity\n");
  printf("\t   1: Mylar (C10-H8-O4)\t\t\t1397 mg/cm3\n");
  printf("\t   2: Polyethylene (C2-H4)\t\t 930 mg/cm3\n");
  printf("\t   3: Polypropylene (C3-H6)\t\t 900 mg/cm3\n");
  printf("\t   4: Kapton (H10-C22-N2-O5)\t\t1430 mg/cm3\n");
  printf("\t   5: Cesium Iodine (CsI)\t\t4510 mg/cm3\n");
  printf("\t   6: Sodium Iodine (NaI)\t\t3670 mg/cm3\n");
  printf("\t   7: Aluminum Oxide (Al2O3)\t\t3980 mg/cm3\n");
  printf("\t   8: Tin-Lead (Sn60/Pb40)\t\t8420 mg/cm3\n");
  printf("\t   9: Nickel (Ni-nat)\t\t\t8902 mg/cm3\n");
  printf("\t  10: Uranium tetrafluoride (UF4)\t6700 mg/cm3\n");
  printf("\t  11: Thorium tetrafluoride (ThF4)\t6300 mg/cm3\n");
  printf("\t  12: Stainless Steel (316L)\t\t7990 mg/cm3\n");
  printf("\n\tGases:\n");
  printf("\t      Gas\t\t\t\tDensity (T=298 K)\n");
  printf("\t 100: Carbon Tetrafluoride (CF4)\tP(Torr) x 4.73538e-3 mg/cm3\n");
  printf("\t 101: Propane (C3-H8)\t\t\tP(Torr) x 2.37276e-3 mg/cm3\n");
  printf("\t 102: Butane (C4-H10)\t\t\tP(Torr) x 3.12752e-3 mg/cm3\n");
  printf("\t 103: Octane (C8-H18)\t\t\tP(Torr) x 6.14657e-3 mg/cm3\n\n");
}

#define ZMAX 118

int zname(char *s) {

  char sym[119][3] = {
    "n" ,"H" ,"He","Li","Be","B" ,"C" ,"N" ,"O" ,"F" ,"Ne","Na","Mg","Al",
    "Si","P" ,"S" ,"Cl","Ar","K" ,"Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co",
    "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y" ,"Zr","Nb",
    "Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I" ,"Xe","Cs",
    "Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
    "Yb","Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi",
    "Po","At","Rn","Fr","Ra","Ac","Th","Pa","U" ,"Np","Pu","Am","Cm","Bk",
    "Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg",
    "Cn","Nh","Fl","Mc","Lv","Ts","Og"};

  for ( int i = 0 ; i <= ZMAX ; i++ ) {
    if ( strcmp(sym[i],s) == 0 ) return i;
  }
  return -1;
}

/*
  Get a sanitized double
*/
bool get_double(double *dest) {
  char s[256];
  scanf("%s",s);
  if ( s == NULL ) {
    return false;
  }
  char *endptr;
  *dest = strtod(s,&endptr);
  if ( s == endptr ) {
    return false;
  }
  // look for trailing text
  while (isspace((unsigned char)*endptr)) {
    endptr++;
  }
  return (*endptr == '\0');
}

void get_ion_data(int *zp, int *ap) {
  char symb[3];
  char cap[256];
  if ( *zp && *ap ) return;
  printf("\n  Please specify data for ion\n");
  printf("\tIon element (e.g. He) : ");
  scanf("%s",symb);
  *zp = zname(&symb[0]);
  if ( *zp == -1 ) {
    fprintf(stderr,"\n\tUnknown element symbol: %s\n\n",symb);
    exit(EXIT_FAILURE);
  }
  printf("\tIon mass number (A)   : ");
  scanf("%s",cap);
  *ap = atoi(cap);
  if ( *ap == 0 ) {
    fprintf(stderr,"\n\tInvalid mass number A: %s\n\n",cap);
    exit(EXIT_FAILURE);
  }
  else if ( *ap < 1 || *ap > 290 ) {
    fprintf(stderr,"\n\tMass number is out of range (0 < A < 290)\n\n");
    exit(EXIT_FAILURE);
  }
}

void get_ion_energy(double *ein) {
  printf("\tGive ion energy (MeV);           %10.2f : ",*ein);
  if ( get_double(ein) ) {
    if ( *ein < 0.0 ) {
      fprintf(stderr,"\n\tInvalid ion energy: %f\n\n",*ein);
      exit(EXIT_FAILURE);
    }
  }
  else {
    fprintf(stderr,"\n\tInvalid entry\n\n");
    exit(EXIT_FAILURE);
  }
}

void get_ion_energy_final(double *eut) {
  printf("\n\tGive final ion energy (MeV);     %10.2f : ",*eut);
  if ( get_double(eut) ) {
    if ( *eut < 0.0 ) {
      fprintf(stderr,"\n\tInvalid final ion energy: %f\n\n",*eut);
      exit(EXIT_FAILURE);
    }
  }
  else {
    fprintf(stderr,"\n\tInvalid entry\n\n");
    exit(EXIT_FAILURE);
  }
}

void get_ion_energy_decrement(double *delen) {
  printf("\n\tGive ion energy decrement (MeV); %10.2f : ",*delen);
  if ( get_double(delen) ) {
    if ( *delen < 0.0 ) {
      fprintf(stderr,"\n\tInvalid ion energy decrement: %f\n\n",*delen);
      exit(EXIT_FAILURE);
    }
  }
  else {
    fprintf(stderr,"\n\tInvalid entry\n\n");
    exit(EXIT_FAILURE);
  }
}

void get_absorber_thickness(double *delr) {
  printf("\n\tGive absorber thickness (mg/cm2);%10.4f : ",*delr);
  if ( get_double(delr) ) {
    if ( *delr < 0.0 ) {
      fprintf(stderr,"\n\tInvalid absorber thickness: %f\n\n",*delr);
      exit(EXIT_FAILURE);
    }
  }
  else {
    fprintf(stderr,"\n\tInvalid entry\n\n");
    exit(EXIT_FAILURE);
  }
}

void make_compound() {
  char symb[3];
  int zel, ael;
  double wel;
  printf("\n  Please specify data for absorber\n");
  printf("\n\tNumber of elements in compound (max %d): ",NELMAX);
  scanf("%d",&nelem);
  if ( nelem < 1 || nelem > NELMAX ) {
    fprintf(stderr,"\n\tIncorrect number of elements\n\n");
    exit(EXIT_FAILURE); 
  }
  printf("\n");
  for ( int i = 0 ; i < nelem ; i++ ) {
    printf("\tElement %1d (e.g. C): ",i+1);
    scanf("%s",symb);
    zel = zname(&symb[0]);
    if ( zel == -1)  {
      fprintf(stderr,"\n\tIllegal element symbol\n\n");
      exit(EXIT_FAILURE);
    }
    printf("\tMass number (A)   : ");
    scanf("%d",&ael);
    if ( ael < 1 || ael > 253 ) {
      fprintf(stderr,"\n\tOut of range (0 < A < 253)\n\n");
      exit(EXIT_FAILURE);
    }
    printf("\tMass weight       : ");
    if ( get_double(&wel) ) {
      if ( wel < 0.0 ) {
	fprintf(stderr,"\n\tInvalid mass weight: %f\n\n",wel);
	exit(EXIT_FAILURE);
      }
    }
    else {
      fprintf(stderr,"\n\tInvalid entry\n\n");
      exit(EXIT_FAILURE);
    }
    absorb[i].z = zel;
    absorb[i].a = ael;
    absorb[i].w = wel;
  }
}

bool isabsorber(int *iabso) {
  if ( *iabso < -1 || *iabso == 0 || *iabso > 103 || (*iabso > 12 && *iabso < 100) ) {
    return false;
  }
  return true;
}

void get_absorber_data(int *iabso, int *zt, int *at) {
  char symb[3];
  char cat[256];
  char answ[256];
  char cabs[256];
  if ( (*zt && *at) || isabsorber(iabso) ) return;
  printf("\n  Please specify data for absorber\n");
  printf("\tIs absorber a (s)ingle element or a (c)ompound [s/c]: ");
  scanf("%s",answ);
  if ( !strcmp(answ,"s") ) {
    printf("\tAbsorber element (e.g. Be) : ");
    scanf("%s",symb);
    *zt = zname(&symb[0]);
    if ( *zt == -1 ) {
      fprintf(stderr,"\n\tIllegal element symbol.\n\n");
      exit(EXIT_FAILURE);
    }
    printf("\tAbsorber mass number (A)   : ");
    scanf("%s",cat);
    *at = atoi(cat);
    if ( *at == 0 ) {
      fprintf(stderr,"\n\tInvalid mass number A: %s\n\n",cat);
      exit(EXIT_FAILURE);
    }
    else if ( *at < 1 || *at > 253 ) {
      fprintf(stderr,"\n\tOut of range (0 < A < 253)\n\n");
      exit(EXIT_FAILURE);
    }
    *iabso = 0;
  }
  else if ( !strcmp(answ,"c") ) {
    list_compounds();
    printf("\tType in the desired absorber: ");
    scanf("%s",cabs);
    *iabso = atoi(cabs);
    if ( *iabso == -1 ) {
      make_compound();
    }
    else if ( !isabsorber(iabso) ) {
      fprintf(stderr,"\n\tInvalid absorber: %s\n\n",cabs);
      exit(EXIT_FAILURE);
    }
  }
  else {
    fprintf(stderr,"\n\tInvalid absorber type: %s\n\n",answ);
    exit(EXIT_FAILURE);
  }
  printf("\n");
}

void check_icorr(int *icorr, int a, double e) {
  if ( *icorr == 0 && e/a > 12.0 ) {
    printf("\n\tE/A > 12 MeV/A. Switching to Hubert-Bimbot-Gauvin correlations\n");
    *icorr = 1;
  }
  if ( *icorr == 1 && e/a <= 2.5 ) {
    printf("\n\tE/A < 2.5 MeV/A. Switching to Northcliffe-Schilling correlations\n");
    *icorr = 0;
  }
}

int main(int argc, char *argv[]) {

  double elin[62] = {0.0100,0.0125,0.0160,0.0200,0.0250,0.0320,0.0400,0.0500,
		     0.0600,0.0700,0.0800,0.0900,0.1000,0.1250,0.1600,0.2000,
		     0.2500,0.3200,0.4000,0.5000,0.6000,0.7000,0.8000,0.9000,
		     1.0000,1.2500,1.5000,2.0000,2.5000,3.2000,4.0000,5.0000,
		     6.0000,7.0000,8.0000,9.0000,10.000,11.000,12.000,15.000,
		     20.000,25.000,30.000,35.000,40.000,45.000,50.000,60.000,
		     70.000,80.000,90.000,100.00,110.00,120.00,150.00,200.00,
		     250.00,300.00,350.00,400.00,450.00,500.00};

  int ilo, ihi;
  double e, tdedxn, tdedxe;
  int n;
  int jj, jjj;
  double err;

  int icorr = 0;
  int zp = 0, ap = 0, zt = 0, at = 0;
  double ein = 0.0, eut = 0.0;
  double delr = 0.0;
  double delen = 0.0;
  int iabso = 0;

  char swt[256];
  int ival;
  char anwr[256];

  // decode command line
  if ( argc > 2 ) {
    disp_header();
    fprintf(stderr,"\nInvalid number of command line arguments.\n\n");
    exit(EXIT_FAILURE);
  }
  else if ( argc == 2 ) {
    if ( !strcmp(argv[1],"-n") || !strcmp(argv[1],"--ns") ) {
      icorr = 0;
    }
    else if ( !strcmp(argv[1],"-h") || !strcmp(argv[1],"--hbg") ) {
      icorr = 1;
    }
    else if ( !strcmp(argv[1],"-l") || !strcmp(argv[1],"--list") ) {
      disp_header();
      printf("\n  List of available compounds:\n");
      list_compounds();
      exit(EXIT_SUCCESS);
    }
    else if ( !strcmp(argv[1],"-v") || !strcmp(argv[1],"--version") ) {
      disp_copy();
      printf("\nrangelib version: %s\n\n",ver);
      exit(EXIT_SUCCESS);
    }
    else if ( !strcmp(argv[1],"--help") ) {
      disp_header();
      disp_help();
      exit(EXIT_SUCCESS);
    }
    else {
      fprintf(stderr,"\nUnknown command line option: %s\n",argv[1]);
      disp_help();
      exit(EXIT_FAILURE);
    }
  }

  disp_header();
  if ( icorr == 0 ) {
    printf("\n  Using the Northcliffe-Schilling correlations (E/A < 12 MeV/A) (default) \n\n");
  }
  if ( icorr == 1 ) {
    printf("\n  Using the Hubert-Bimbot-Gauvin correlations (2.5 < E/A < 100 MeV/A)\n\n");
  }
  printf("  You can choose between the following options:\n\n");
  printf("\t1: Store a table of -dE/dx values\n");
  printf("\t2: Store a table of range values\n");
  printf("\t3: Display the energy after passage of a foil\n");
  printf("\t4: Display foil thickness for a given energy decrement\n");
  printf("\t5: Display foil thickness for a given final energy\n");
  printf("\t6: Display range for a given initial energy\n");
  printf("\t7: Display initial energy from a given range\n");
  printf("\t8: Display initial energy for a given final energy\n");
  printf("\t9: Exit\n");

  printf("\n  Type in the desired option: ");
  scanf("%s",swt);
  ival = atoi(swt);
  if ( ival < 1 || ival > 9 ) {
    fprintf(stderr,"\n  Invalid option: %s\n\n",swt);
    exit(EXIT_FAILURE);
  }

  while ( 1 ) {
    switch(ival) {
    case 1:
      printf("\tGive a file name for the resulting data: ");
      char filnamn[256];
      scanf("%s",filnamn);
      FILE *fp = fopen(filnamn,"w");
      get_ion_data(&zp,&ap);
      get_absorber_data(&iabso,&zt,&at);
      fprintf(fp,"   E/A         E      -dE/dx(el)  -dE/dx(nu) -dE/dx(tot)\n");
      fprintf(fp," (MeV/A)     (MeV)   (MeVcm2/mg) (MeVcm2/mg) (MeVcm2/mg)\n");
      if ( icorr == 0 ) {
	ilo = 0;
	ihi = 39;
      }
      else {
	ilo = 28;
	ihi = 62;
      }
      for ( int i = ilo ; i < ihi ; i++ ) {
	e = elin[i];
	dedxtab(icorr,zp,ap,iabso,zt,at,e,&tdedxe,&tdedxn);
	fprintf(fp,"%8.4f %10.4f %11.4f %11.4f %11.4f\n",
		e,e*ap,tdedxe,tdedxn,tdedxn+tdedxe);
      }
      fclose(fp);
      break;
    case 2:
      printf("\tGive a file name for the resulting data: ");
      scanf("%s",filnamn);
      fp = fopen(filnamn,"w");
      get_ion_data(&zp,&ap);
      get_absorber_data(&iabso,&zt,&at);
      double *em = malloc(NMAX*sizeof(double));
      double *r = malloc(NMAX*sizeof(double));
      rangetab(icorr,zp,ap,iabso,zt,at,&em[0],&r[0],&n);
      jjj = n-3;
      fprintf(fp,"   E/A         E           R\n");
      fprintf(fp," (MeV/A)     (MeV)      (mg/cm2)\n");
      if ( icorr == 0 ) {
	ilo = 0;
	ihi = 39;
      }
      else {
	ilo = 28;
	ihi = 62;
      }
      for ( int i = ilo ; i < ihi ; i++ ) {
	e = log10(elin[i]);
	jj = nr_locate(em,n,e);
	if ( jj > jjj ) jj = jjj;
	fprintf(fp,"%8.4f %10.4f %11.4f\n",
		elin[i],elin[i]*ap,nr_polint(&em[jj],&r[jj],3,e,&err));
      }
      // free allocated memory
      free(em);
      free(r);
      fclose(fp);
      break;
    case 3:
      get_ion_data(&zp,&ap);
      get_absorber_data(&iabso,&zt,&at);
      get_ion_energy(&ein);
      check_icorr(&icorr,ap,ein);
      get_absorber_thickness(&delr);
      eut = passage(icorr,zp,ap,iabso,zt,at,ein,delr,&err);
      if ( eut == 0.0 ) {
	printf("\n\tAbsorber will stop all ions\n\n");
      }
      else {
	printf("\n\tInitial ion energy : %10.4f MeV\n",ein);
	printf("\tFinal ion energy   : %10.4f MeV\n",eut);
	printf("\tEnergy loss in foil: %10.4f MeV\n\n",ein-eut);
      }
      break;
    case 4:
      get_ion_data(&zp,&ap);
      get_absorber_data(&iabso,&zt,&at);
      get_ion_energy(&ein);
      check_icorr(&icorr,ap,ein);
      get_ion_energy_decrement(&delen);
      if ( delen > ein ) {
	printf("\n\tCalculation is not possible. Energy decrement > initial energy\n\n");
      }
      else {
	delr = thickn(icorr,zp,ap,iabso,zt,at,ein,delen);
	printf("\n\tInitial ion energy : %10.4f MeV\n",ein);
	printf("\tFinal ion energy   : %10.4f MeV\n",ein-delen);
	printf("\tFoil thickness is  : %10.4f mg/cm2\n\n",delr);
      }
      break;
    case 5:
      get_ion_data(&zp,&ap);
      get_absorber_data(&iabso,&zt,&at);
      get_ion_energy(&ein);
      check_icorr(&icorr,ap,ein);
      get_ion_energy_final(&eut);
      if ( eut > ein ) {
	printf("\n\tCalculation is not possible: Final energy > initial energy\n\n");
      }
      else {
	delr = thickn(icorr,zp,ap,iabso,zt,at,ein,ein-eut);
	printf("\n\tInitial ion energy : %10.4f MeV\n",ein);
	printf("\tFinal ion energy   : %10.4f MeV\n",eut);
	printf("\tFoil thickness is  : %10.4f mg/cm2\n\n",delr);
      }
      break;
    case 6:
      get_ion_data(&zp,&ap);
      get_absorber_data(&iabso,&zt,&at);
      get_ion_energy(&ein);
      check_icorr(&icorr,ap,ein);
      delr = rangestop(icorr,zp,ap,iabso,zt,at,ein);
      printf("\n\tInitial ion energy : %10.4f MeV\n",ein);
      printf("\tThe range is       : %10.4f mg/cm2\n\n",delr);
      break;
    case 7:
      get_ion_data(&zp,&ap);
      get_absorber_data(&iabso,&zt,&at);
      get_absorber_thickness(&delr);
      ein = rangein(icorr,zp,ap,iabso,zt,at,delr);
      int isav = icorr;
      check_icorr(&icorr,ap,ein);
      if ( icorr != isav ) {
	ein = rangein(icorr,zp,ap,iabso,zt,at,delr);
      }
      printf("\n\tMaximum ion energy stopped by this foil: %10.4f MeV\n\n",ein);
      break;
    case 8:
      get_ion_data(&zp,&ap);
      get_absorber_data(&iabso,&zt,&at);
      get_ion_energy_final(&eut);
      get_absorber_thickness(&delr);
      ein = egassap(icorr,zp,ap,iabso,zt,at,delr,eut,&err);
      printf("\n\tInitial ion energy : %10.4f MeV\n",ein);
      printf("\tFinal ion energy   : %10.4f MeV\n",eut);
      printf("\tEnergy loss in foil: %10.4f MeV\n\n",ein-eut);
      break;
    case 9:
      exit(EXIT_SUCCESS);
      break;
    }
    printf("  Another calculation? (y/n): ");
    scanf("%s",anwr);
    if ( !strcmp(anwr,"n") ) {
      printf("\n");
      break;
    }
    else if ( strcmp(anwr,"y") ) {
      fprintf(stderr,"\n  Invalid input: %s\n\n",anwr);
      exit(EXIT_FAILURE);
    }
    else {
      printf("\n");
    }
  }
  return EXIT_SUCCESS;
}
