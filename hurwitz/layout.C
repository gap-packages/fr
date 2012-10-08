// solve Hurwitz problem by laying out triangulation

#undef PRINT_DOGLEG
#undef PRINT_U
#undef PRINT_XY
#define PRINT_STATS

#include "config.h"
#include <iostream>
#include <queue>
#include <map>
using namespace std;
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <cholmod.h>
extern "C" {
#include "dogleg.h"
}

#ifndef HAVE_SINCOS
void inline sincos(double x, double *s, double *c)
{
  *s = sin(x);
  *c = cos(x);
}
#endif

cholmod_common cholmod;

#define MAXVERTICES 100000
#define MAXFACES (2*(MAXVERTICES-2))
#define maxiter 1000

int numvertices, innervertices, ivertex, numfaces, jac_nnz;
int vertex[MAXFACES][3]; // vertices of face
int face[MAXFACES][3]; // face[f][i] is face opposite vertex[f][i]
double loglength[MAXFACES][3]; // log-lengths of edges opposite to vertices
bool boundary[MAXVERTICES];
double u[MAXVERTICES]; // will hold the scaling at vertices
double x[MAXVERTICES], y[MAXVERTICES]; // laid-out positions on the plane

// inner2v enumerates the inner vertices, v2inner gives the inner# of a vertex
int v2inner[MAXVERTICES], inner2v[MAXVERTICES];

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

int numbrokentriangles; // statistics

void triangle_angles (double loglen[3], double angles[3], double cot[3])
{
  double len[3];
  for (int i = 0; i < 3; i++) len[i] = exp(loglen[i]);
  double s0 = len[0] + len[1] + len[2];
  double s[3];
  for (int i = 0; i < 3; i++) s[i] = s0 - 2*len[i];

  if (s[0] <= 0. || s[1] <= 0. || s[2] <= 0.) {
    numbrokentriangles++;
    for (int i = 0; i < 3; i++) angles[i] = M_PI * (s[i] <= 0.), cot[i] = 0.;
    return;
  }
  for (int i = 0; i < 3; i++)
    angles[i] = 2. * atan(sqrt(s[(1+i)%3] * s[(2+i)%3] / (s[i] * s0)));
  double p = 0.5 / sqrt(s[0] * s[1] * s[2] * s0);
  for (int i = 0; i < 3; i++)
    cot[i] = p * (s[i] * s0 - s[(1+i)%3] * s[(2+i)%3]);
}

void add_triplet (cholmod_triplet *matrix, int i, int j, double x)
{
  ((int *)matrix->i)[matrix->nnz] = i;
  ((int *)matrix->j)[matrix->nnz] = j;
  ((double *)matrix->x)[matrix->nnz] = x;
  matrix->nnz++;
}

void layout_fdf (const double x[], double y[], cholmod_sparse *jacobian, void* cookie __attribute__ ((unused)))
{
  for (int i = 0; i < innervertices; i++) // fill in u values
    u[inner2v[i]] = x[i];

  numbrokentriangles = 0; // for statistics

  cholmod_triplet *j_triplet = cholmod_allocate_triplet (innervertices, innervertices, jac_nnz, 0, CHOLMOD_REAL, &cholmod); // in fact, is symmetric

  double value = 0., gradient[numvertices];

  for (int i = 0; i < numvertices; i++)
    value += 2.*M_PI*u[i];

  for (int i = 0; i < innervertices; i++)
    y[i] = 2.*M_PI;

  for (int i = 0; i < numfaces; i++) {
    double newloglength[3], angle[3], cot[3];

    for (int j = 0; j < 3; j++)
      newloglength[j] = loglength[i][j] + u[vertex[i][(j+1)%3]] + u[vertex[i][(j+2)%3]];
    triangle_angles(newloglength, angle, cot);

    for (int j = 0; j < 3; j++) {
      value -= M_PI*u[vertex[i][j]];
      value += angle[j]*newloglength[j];
      value += 0.5 * clausen(2.*angle[j]);

      int vertexij = v2inner[vertex[i][j]];
      if (vertexij == -1)
	continue;

      y[vertexij] -= angle[j];

      add_triplet (j_triplet, vertexij, vertexij, cot[(j+1)%3] + cot[(j+2)%3]);
      for (int k = 0; k < 3; k++) {
	if (j == k) continue;
	int vertexik = v2inner[vertex[i][k]];
	if (vertexik == -1)
	  continue;
	add_triplet (j_triplet, vertexij, vertexik, -cot[3-j-k]);
      }
    }
  }

  cholmod_sparse *jac = cholmod_triplet_to_sparse (j_triplet, jac_nnz, &cholmod);
  memcpy (jacobian->p, jac->p, (jac->ncol+1)*sizeof(int));
  memcpy (jacobian->i, jac->i, jac->nzmax*sizeof(int));
  memcpy (jacobian->x, jac->x, jac->nzmax*sizeof(double));

  cholmod_free_sparse (&jac, &cholmod);
  cholmod_free_triplet (&j_triplet, &cholmod);

#ifdef PRINT_STATS
  cerr.precision(16);
  cerr << numbrokentriangles << " broken triangles, value = " << value << endl;
#endif
}

void computelengths (void)
{
  innervertices = 0;
  for (int i = 0; i < numvertices; i++)
    if (boundary[i])
      v2inner[i] = -1;
    else {
      v2inner[i] = innervertices;
      inner2v[innervertices++] = i;
    }
  jac_nnz = 9*numfaces;

  cerr << "infty = " << ivertex << ", " << numvertices-innervertices << " boundary vertices, " << innervertices << " interior vertices." << endl;

  double inneru[innervertices];
  for (int i = 0; i < innervertices; i++)
    inneru[i] = 0.;

#ifdef PRINT_DOGLEG
  dogleg_setDebug(1);
#else
  dogleg_setDebug(0);
#endif
  dogleg_setMaxIterations (1000);
  dogleg_setThresholds (1.e-12,0.,-1.);

#ifdef CHECK_GRADIENT
  fprintf(stderr, "have %d variables\n", innervertices);
  for(int i=0; i<innervertices; i++)
  {
    fprintf(stderr, "checking gradients for variable %d\n", i);
    dogleg_testGradient(i, u, innervertices, innervertices, jac_nnz, layout_fdf, NULL);
  }
#endif
 
  double optimum = dogleg_optimize(inneru, innervertices, innervertices, jac_nnz, layout_fdf, NULL, NULL);

#ifdef PRINT_DOGLEG
  cerr << "Achieved norm " << optimum << endl;
#endif

  for (int i = 0; i < innervertices; i++)
    u[inner2v[i]] = inneru[i];
}

void makefaces (void)
{
  typedef pair<int,int> edge;

  map<edge,int> edge2triangle;

  for (int i = 0; i < numfaces; i++)
    for (int j = 0; j < 3; j++)
      edge2triangle.insert (pair<edge,int>(edge(vertex[i][j],vertex[i][(j+1)%3]),i));

  for (int i = 0; i < numfaces; i++)
    for (int j = 0; j < 3; j++)
      face[i][j] = edge2triangle[edge(vertex[i][(j+2)%3],vertex[i][(j+1)%3])];
}

void layoutfaces (void)
{
  queue<int> q;
  bool vtag[numvertices], ftag[numfaces];
  double length[numfaces][3];

  for (int i = 0; i < numvertices; i++)
    vtag[i] = false;
  for (int i = 0; i < numfaces; i++)
    ftag[i] = false;
  vtag[ivertex] = true;

  for (int i = 0; i < numfaces; i++)
    for (int j = 0; j < 3; j++)
      length[i][j] = exp(loglength[i][j] + u[vertex[i][(j+1)%3]] + u[vertex[i][(j+2)%3]]);

  int rootface;
  for (rootface = 0;; rootface++) {
    int j;
    for (j = 0; j < 3; j++)
      if (vertex[rootface][j] == ivertex) break;
    if (j == 3) break;
  }
  q.push (rootface);

  while (!q.empty()) {
    int f = q.front();
    q.pop();
    if (ftag[f])
      continue;
    ftag[f] = true;
    for (int j = 0; j < 3; j++) {
      q.push(face[f][j]);

      int v0 = vertex[f][j];
      if (vtag[v0])
	continue;

      // v0 has not been laid. Lay it now.

      vtag[v0] = true;
      int v1 = vertex[f][(j+1)%3], v2 = vertex[f][(j+2)%3];

      if (!vtag[v1] && !vtag[v2]) { // very first vertex
	x[v0] = y[v0] = 0.;
	continue;
      }

      if (!vtag[v1] || !vtag[v2]) { // very second vertex
	x[v0] = length[f][(j+1)%3];
	y[v0] = 0.;
	continue;
      }
      // now v1 and v2 have already been laid out. compute pos. of v0

      double e0 = length[f][j],
	e1 = length[f][(j+1)%3],
	e2 = length[f][(j+2)%3];
      double angle1 = (e0*e0 + e2*e2 - e1*e1) / (2.*e0*e2);
      if (angle1 >= 1.0) // safety, for roundoffs
	angle1 = 0.0;
      else if (angle1 <= -1.0)
	angle1 = M_PI;
      else
	angle1 = acos(angle1);

      double slope12 = atan2(y[v2]-y[v1],x[v2]-x[v1]);
      double dx, dy;
      sincos (slope12+angle1, &dy, &dx);
      x[v0] = x[v1] + e2*dx;
      y[v0] = y[v1] + e2*dy;
    }
  }
#if 0 // sanity check
  for (int i = 0; i < numfaces; i++)
    for (int j = 0; j < 3; j++) {
      int v0 = vertex[i][(j+1)%3], v1 = vertex[i][(j+2)%3];
      cerr << "face " << i << " edge " << j << ": length " << hypot(x[v0]-x[v1],y[v0]-y[v1]) << " (expected " << length[i][j] << ")\n";
    }
#endif

  for (int i = 0; i < numvertices; i++)
    if (!vtag[i]) {
      cerr << "Didn't tag vertex " << i << ". Repent.\n";
      exit(-1);
    }
}

int main(int argc, char *argv[])
{
  cholmod_start (&cholmod);

  for (;;) {
    char s[100];
    if (scanf("%s", s) != 1)
      strcpy(s,"EOF");
    if (!strcmp(s,"END"))
      break;
    else if (!strcmp(s,"VERTICES")) {
      scanf("%d", &numvertices);
      scanf("%d", &ivertex);
      ivertex--; // C array indexing
      if (numvertices > MAXVERTICES) {
	cerr << "numvertices = " << numvertices << " > MAXVERTICES = " << MAXVERTICES << ". Repent." << endl;
	return -1;
      }
    } else if (!strcmp(s,"FACES")) {
      scanf ("%d", &numfaces);
      for (int i = 0; i < numfaces; i++) {
	for (int j = 0; j < 3; j++) {
	  scanf("%d", vertex[i]+j);
	  vertex[i][j]--; // C arrays start at 0
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
      if (vertex[i][j] == ivertex) {
	boundary[vertex[i][(j+1)%3]] = true;
	u[vertex[i][(j+1)%3]] = -loglength[i][(j+2)%3];
      }

  computelengths();

#ifdef PRINT_U
  cerr << "u vector:\n--------\n";
  for (int i = 0; i < numvertices; i++)
    cerr << u[i] << endl;
#endif

  cholmod_finish (&cholmod);

  makefaces ();
  layoutfaces ();

#ifdef PRINT_XY
  cerr << "points:\n------\n";
  for (int i = 0; i < numvertices; i++)
    cerr << "(" << x[i] << "," << y[i] << ")" << endl;
#endif

  // center the points
  double sx = 0., sy = 0., s = 0.;
  for (int i = 0; i < numvertices; i++)
    sx += x[i], sy += y[i];
  sx /= numvertices; sy /= numvertices;
  for (int i = 0; i < numvertices; i++)
    x[i] -= sx, y[i] -= sy;

  // spread them out
  for (int i = 0; i < numvertices; i++)
    s += x[i]*x[i]+y[i]*y[i];
  s = sqrt(s / numvertices);

  for (int i = 0; i < numvertices; i++)
    x[i] /= s, y[i] /= s;

  if (s != s) // we got NaN, because one of the values is invalid
    return -1;
    
  // print stereographic projection
  cout.precision(20);
  cout << "[";
  for (int i = 0; i < numvertices; i++) {
    double n = x[i]*x[i] + y[i]*y[i] + 1.;

    if (i == ivertex) {
      x[i] = y[i] = 0.;
      n = 1.e300; // just as good as infinity, if we don't have 300 sig.digits
    }

    cout << "[" << 2.*x[i]/n << "," << 2.*y[i]/n << "," << (n-2.)/n << "]," << endl;
  } 
  cout << "fail];" << endl;

  return 0;
}
