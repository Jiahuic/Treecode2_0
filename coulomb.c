/*
 * coulomb.c: main driver
 * This program computes the screened Coulomb potential with treecode method
 *         exp(-K|x-y|)
 * \phi = --------------
 *            |x-y|
 * usage:
 *   coulomb [options] panelfile [options]
 *
 * Copyright Jiahui Chen
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "treeGlobal.h"
#include "treecode.h"
#include <time.h>

#define MEG 1048576

/* routines used by the main routine */
double *loadPanel(char *panelfile, char *density, int *numSing);
void setupCOEF(ssystem *sys);
void createTree(cube *cb);
void compMomAll(cube *cb, int ifirst);
void compTree(cube *cb, double *pot);
void cubeNumber(ssystem *sys);
void directSum(ssystem *sys, double *dpot);
void printError(ssystem *sys, double *pot, double *dpot);
void printMemory(ssystem *sys);
void freeCOEF(ssystem *sys, double *pot, double *dpot);

int main(int nargs, char *argv[]){
  char panelfile[80], density[80];
  int i, j, k, n, nPnts;
  double *pot, *dpot, abserr, relerr;
  ssystem *sys;

  /* time variables */
  clock_t start, diff;
  int msec;

  /* set up system */
  CALLOC(sys, 1, ssystem, ON, AMISC);
  sprintf(density, "3");
  sprintf(panelfile, "1a63");
  sys->kappa = 1.0;
  sys->theta = 0.8;
  sys->order = 3;
  sys->maxparCube = 500;

  for ( i=1; i<nargs; i++ )
    if ( argv[i][0] == '-' )
      switch ( argv[i][1] ) {
      case 'p': sys->order = atoi( argv[i]+3 );
        break;
      case 'r': sys->theta = atof( argv[i]+3 );
        break;
      case 'd': strcpy(density,argv[i]+3);
        break;
      case 'n': sys->maxparCube = atoi( argv[i]+3 );
        break;
      }
    else
      strcpy(panelfile,argv[i]);

  sys->positions = loadPanel(panelfile, density, &nPnts);
  sys->nPnts = nPnts;
  CALLOC(sys->topCube, 1, cube, ON, ACUBES);
  CALLOC(sys->charges, nPnts, double, ON, AXYZQ);
  for ( i=0; i<nPnts; i++ ) sys->charges[i] = 1.0;

  CALLOC(pot, nPnts, double, ON, AMISC);
  CALLOC(dpot, nPnts, double, ON, AMISC);


  /* setup treecode coefficients */
  setupCOEF(sys);

  /* build tree */
  start = clock();
  createTree(sys->topCube);
  diff = clock() - start;
  msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Treecode create Tree time %d sec %d msec\n", msec/1000, msec%1000);

  start = clock();
  compMomAll(sys->topCube, 1);
  diff = clock() - start;
  msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Moments computation time %d sec %d msec\n", msec/1000, msec%1000);

  /* compute treecode, output as "pot" */
  start = clock();
  compTree(sys->topCube, pot);
  diff = clock() - start;
  msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Treecode compute Tree time %d sec %d msec\n", msec/1000, msec%1000);

  /* compute direct sum */
  start = clock();
  directSum(sys, dpot);
  diff = clock() - start;
  msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("Direct sum time %d sec %d msec\n", msec/1000, msec%1000);

  /* error estimate, then print memory usage */
  printf("  \n");
  printf("Run parameteres:  \n");
  printf("                   numpars    = %d\n",sys->nPnts);
  printf("                   kappa      = %f\n",sys->kappa);
  printf("                   theta      = %f\n",sys->theta);
  printf("                   order      = %d\n",sys->order);
  printf("                   maxparnode = %d\n",sys->maxparCube);
  printf("                   max level  = %d\n",sys->maxLev);
  printError(sys, pot, dpot);
  printMemory(sys);

  freeCOEF(sys, pot, dpot);
  free(sys->positions);
  free(sys->charges);
  free(pot);
  free(dpot);

  free(sys->topCube);
  free(sys);

}

void directSum(ssystem *sys, double *dpot) {
  int i, j, k, nPnts=sys->nPnts;
  double *source, *target, dist, dx, dy, dz, energy, temp;
  double kappa=sys->kappa;

  for ( i=0; i<nPnts; i++ ){
    energy = 0.0;
    target = &sys->positions[3*i];
    for ( j=i+1; j<nPnts; j++ ){
      source = &sys->positions[3*j];

      dx = source[0]-target[0];
      dy = source[1]-target[1];
      dz = source[2]-target[2];
      dist = sqrt(SQR(dx)+SQR(dy)+SQR(dz));
      //temp = exp(-kappa*dist)/dist*fourPiI;
      temp = exp(-kappa*dist)/dist;

      energy += sys->charges[j]*temp;
      dpot[j] += sys->charges[i]*temp;
    }
    dpot[i] = sys->charges[i]*(dpot[i]+energy);
  }
  dpot[nPnts] = sys->charges[nPnts]*dpot[nPnts];
} /* directSum */


void printError(ssystem *sys, double *pot, double *dpot) {
  int i, nPnts=sys->nPnts;
  double relerr, abserr;

  relerr = 0.0; abserr = 0.0;
  for ( i=0; i<nPnts; i++ ){
    relerr += SQR(pot[i]-dpot[i]);
    abserr += SQR(dpot[i]);
  }
  relerr = sqrt(relerr/abserr);

  printf("the relative err: %.15f\n",relerr);
} /* printError */

void printMemory(ssystem *sys) {
  double ratio;
  long total;
  cube *cb;

  printf("\nNumber points=%d", sys->nPnts);
  printf("\nMemory: Total                 %lg MB\n", ((double)memcount)/MEG);
  printf("        Positions/Charges     %lg MB\n", ((double)memXYZQ)/MEG);
  printf("        Tree                  %lg MB\n", ((double)memCUBES)/MEG);
  printf("        System                %lg MB\n", ((double)memMISC)/MEG);

}

void cubeSelfCount(cube *cb, int counter){
  int i;
  cube *kid;

  for ( i=0; i<cb->nKids; i++ ) {
    kid = cb->kids[i];
    cubeSelfCount(kid, counter);
    counter++;
  }
}

/* count cube number for each level */
void cubeNumber(ssystem *sys) {
  int i, lev, cubeN, tcubeN=0;
  cube *cb=sys->topCube;

  cubeSelfCount(cb, tcubeN);
  printf("Total cube number %d\n",tcubeN);
}
