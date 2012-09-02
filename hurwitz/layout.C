// solve Hurwitz problem by laying out triangulation

#include <iostream>
using namespace std;
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "levmar.h"

#define MAXVERTICES 10000
#define MAXFACES (2*(MAXVERTICES-2))
#define maxiter 1000

int numvertices, innervertices, ivertex, numfaces;
int face[MAXFACES][3]; // vertices of face
double loglength[MAXFACES][3]; // log-lengths of edges opposite to vertices
bool boundary[MAXVERTICES];

// inner2v enumerates the inner vertices, v2inner gives the inner# of a vertex
int v2inner[MAXVERTICES], inner2v[MAXVERTICES];

double u[MAXVERTICES]; // will hold the scaling at vertices

double clausen (double x)
// Clausen's integral
{
  x = fmod(x + M_PI, 2.*M_PI) - M_PI;

  if (x == 0.)
    return 0.;

  if (fabs(x) <= 2.0944) { // small
    double xx = x*x;

    return ((((((((((((2.3257441143020875e-22 * xx
		       + 1.0887357368300848e-20) * xx
		      + 5.178258806090624e-19) * xx
		     + 2.5105444608999545e-17) * xx
		    + 1.2462059912950672e-15) * xx
		   + 6.372636443183181e-14) * xx
		  + 3.387301370953521e-12) * xx
		 + 1.8978869988971e-10) * xx
		+ 1.1482216343327455e-8) * xx
	       + 7.873519778281683e-7) * xx
	      + 0.00006944444444444444) * xx
	     + 0.013888888888888888) * xx
	    - log(fabs(x)) + 1.0) * x;
  } else { // big
   x -= M_PI*((x > 0.) - (x < 0.));
   double xx = x*x;
   return ((((((((((((3.901950904063069e-15 * xx
		      + 4.566487567193635e-14) * xx
		     + 5.429792727596476e-13) * xx
		    + 6.5812165661369675e-12) * xx
		   + 8.167010963952222e-11) * xx
		  + 1.0440290284867003e-9) * xx
		 + 1.3870999114054669e-8) * xx
		+ 1.941538399871733e-7) * xx
	       + 2.927965167548501e-6) * xx
	      + 0.0000496031746031746) * xx
	     + 0.0010416666666666667) * xx
	    + 0.041666666666666664) * xx
	   + -0.693147180559945) * x;
  }
}

int numbrokentriangs; // statistics

void triangle_angles (double loglen[3], double angles[3], double cot[3])
{
  double len[3];
  for (int i = 0; i < 3; i++) len[i] = exp(loglen[i]);
  double s0 = len[0] + len[1] + len[2];
  double s[3];
  for (int i = 0; i < 3; i++) s[i] = s0 - 2*len[i];

  if (s[0] <= 0. || s[1] <= 0. || s[2] <= 0.) {
    numbrokentriangs++;
    for (int i = 0; i < 3; i++) angles[i] = M_PI * (s[i] <= 0.), cot[i] = 0.;
    return;
  }
  for (int i = 0; i < 3; i++)
    angles[i] = 2. * atan(sqrt(s[(1+i)%3] * s[(2+i)%3] / (s[i] * s0)));
  double p = 0.5 / sqrt(s[0] * s[1] * s[2] * s0);
  for (int i = 0; i < 3; i++)
    cot[i] = p * (s[i] * s0 - s[(1+i)%3] * s[(2+i)%3]);
}

double lastfunctionvalue, thisfunctionvalue;

void layout_fdf (double x[], double y[], double z[])
{
  for (int i = 0; i < innervertices; i++) // fill in u values
    u[inner2v[i]] = x[i];

  numbrokentriangs = 0;

  double value = 0., gradient[numvertices];

  for (int i = 0; i < numvertices; i++) {
    value += 2.*M_PI*u[i];
    gradient[i] = 2.*M_PI;
  }

  for (int i = 0; i < numfaces; i++) {
    double newloglength[3], angle[3], cot[3];

    for (int j = 0; j < 3; j++)
      newloglength[j] = loglength[i][j] + u[face[i][(j+1)%3]] + u[face[i][(j+2)%3]];
    triangle_angles(newloglength, angle, cot);

    for (int j = 0; j < 3; j++) {
      value -= M_PI*u[face[i][j]];
      value += angle[j]*newloglength[j];
      value += 0.5 * clausen(2.*angle[j]);
    }

    lastfunctionvalue = thisfunctionvalue;
    thisfunctionvalue = value;

    for (int j = 0; j < 3; j++)
      gradient[face[i][j]] -= angle[j];

#ifdef WITH_HESSIAN
    for (int j = 0; j < 3; j++) {
      hessian[face[i][j]][face[i][j]] = cot[(j+1)%3] + cot[(j+2)%3];
      for (int k = 0; k < 3; k++) {
	if (j == k) continue;
	hessian[face[i][j]][face[i][k]] = -cot[3-j-k];
      }
    }

    for (int i = 0; i < innervertices; i++)
      for (int j = 0; j < innervertices; j++)
	hess[i][j] = hessian[inner2v[i]][inner2v[j]];
#endif
  }

  if (y)
    y[0] = value;

  if (z)
    for (int i = 0; i < innervertices; i++)
      z[i] = gradient[inner2v[i]];
}

void layout_f (double *p, double *hx, int m, int n, void *data)
{
  layout_fdf (p, hx, NULL);
}

void layout_df (double *p, double *j, int m, int n, void *data)
{
  layout_fdf (p, NULL, j);
}

void printinfo (int iter, int call_iter, double x[], double *f, double *g,  double* gnorm)
{
  cout << iter << ": " << call_iter << " " << *f << " " << *gnorm  << " broken: " << numbrokentriangs << "... \r";
}

int main(int argc, char *argv[])
{
  for (;;) {
    char s[100];
    if (scanf("%s", s) != 1)
      break;
    if (!strcmp(s,"VERTICES")) {
      scanf("%d", &numvertices);
      scanf("%d", &ivertex);
      ivertex--; // C array indexing
    } else if (!strcmp(s,"FACES")) {
      scanf ("%d", &numfaces);
      for (int i = 0; i < numfaces; i++) {
	for (int j = 0; j < 3; j++) {
	  scanf("%d", face[i]+j);
	  face[i][j]--; // C arrays start at 0
	}
	for (int j = 0; j < 3; j++) {
	  double d;
	  scanf("%lf", &d);
	  loglength[i][j] = log(d);
	}
      }
    } else {
      cerr << "Unknown input " << s << ". Repent." << endl;
      return -1;
    }
  }

  if (numfaces != 2*(numvertices-2)) {
    cerr << "Error: numvertices = " << numvertices << ", numfaces = " << numfaces << ". Repent." << endl;
    return -1;
  }

  for (int i = 0; i < numvertices; i++)
    boundary[i] = false;
  boundary[ivertex] = true;

  for (int i = 0; i < numfaces; i++)
    for (int j = 0; j < 3; j++)
      if (face[i][j] == ivertex) {
	boundary[face[i][(j+1)%3]] = true;
	u[face[i][(j+1)%3]] = -loglength[i][(j+2)%3];
      }

  innervertices = 0;
  for (int i = 0; i < numvertices; i++)
    if (boundary[i])
      v2inner[i] = -1;
    else {
      v2inner[i] = innervertices;
      inner2v[innervertices++] = i;
    }

  cerr << "infty = " << ivertex << ", " << numvertices-innervertices << " boundary vertices, " << innervertices << " interior vertices." << endl;

  double inneru[innervertices];
  for (int i = 0; i < innervertices; i++)
    inneru[i] = 0.;
  double min;
  double info[LM_INFO_SZ];

  int iter = dlevmar_der (layout_f, layout_df, inneru, &min, innervertices, 1, 100, NULL, info, NULL, NULL, NULL);

  if (iter < 0) {
    cerr << "I freaked out." << endl;
    return -1;
  }

  cerr << min << endl;

  switch ((int) info[6]) {
  case 1:
    cerr << "stopped by small gradient J^T e" << endl;
    break;
  case 2:
    cerr << "stopped by small Dp" << endl;
    break;
  case 3:
    cerr << "stopped by itmax" << endl;
    break;
  case 4:
    cerr << "singular matrix. Restart from current p with increased μ" << endl;
    break;
  case 5:
    cerr << "no further error reduction is possible. Restart with increased μ" << endl;
    break;
  case 6:
    cerr << "success after " << iter << " iterations, " << info[7] << " function evals, " << info[8] << " jacobians, " << info[9] << " linear systems." << endl;
    break;
  case 7:
    cerr << "stopped by invalid (i.e. NaN or Inf) function values. Repent." << endl;
    return -1;
  }

  cerr << "|e|₂ = " << info[1] << ", |J⁺e|ₒₒ = " << info[2] << ", |Dp|₂ = " << info[3] << ", μ/max[J⁺J]ₖₖ ] = " << info[4] << endl << endl;

  for (int i = 0; i < innervertices; i++)
    u[inner2v[i]] = inneru[i];
}
