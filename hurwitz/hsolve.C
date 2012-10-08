// solve Hurwitz's problem, by Laurent Bartholdi, 20120827

// requires levmar library for Levenberg-Marquardt,
// see http:// www.ics.forth.gr/~lourakis/levmar/

#include <iostream>
#include <complex>
using namespace std;
#include <string.h>
#include "levmar.h"

typedef complex<double> cdouble;

#define MAXDEGREE 100

const int maxiter = 2000;

// rational maps are normalized so that 0 -> 0, 1 -> 1, infty -> infty;
// so they have the form A\prod_{i=0}^numzero (z-p_i)^e_i for some e_i in Z\{0},
// and we force c_numzero = 0 and e_numzero > 0 and A such that
// prod_{i=0}^numzero (1-p_i)^e_i = 1.
// they also have critical points c_0,...,c_{numcv-1},
// mapping to critical values v_0,...,v_{numcv-1} respectively.
// critical value c_i has multiplicity d_i.

struct hurwitzdata {
  int degree, // degree of rational map
    numzero, // number of zeros and poles, excepting 0,infty
    numcv; // number of critical values, excepting poles,zeros and 1

  int e[2*MAXDEGREE-1], // orders of zeros and poles, + for zero, - for pole
                        // e[numzero] is the order of the zero at 0
    d[2*MAXDEGREE-1]; // orders of critical points.
                      // d[numcv] is the order of 1
  cdouble v[2*MAXDEGREE-2]; // desired critical values
};

// computes z^i for i integer
cdouble ipow (cdouble z, int i)
{
  cdouble v = 1.0;

  if (i < 0) i = -i, z = 1.0/z;
  for (;;) {
    if (i & 1) v *= z;
    if (i >>= 1) z *= z; else break;
  }
  return v;
}

// not really a distance, it's complex-valued; but a good measure of
// how points differ on the sphere
cdouble p1distance (cdouble u, cdouble v)
{
  return (u-v) / sqrt((1.0+norm(u)) / (1.0+norm(v)));
}

void solvehd_f (cdouble x[], cdouble y[], int m, int n, const hurwitzdata &data)
{
  cdouble p[data.numzero+1], c[data.numcv+1];
  int eq = 0;

  // the equations we take are:
  // at each critical point c[i], (f'/f)^(k)=0 for k=0,...,data.d[i]-1
  // at each critical point c[i], p1distance(f(c[i])/f(1),data.v[i])=0

  for (int i = 0; i < data.numzero; i++)
    p[i] = x[i];
  p[data.numzero] = 0.0; // last zero is always 0

  for (int i = 0; i < data.numcv; i++)
    c[i] = x[i+data.numzero];
  c[data.numcv] = 1.0; // last critical point is always 1

  cdouble A = 1.0;
  for (int i = 0; i < data.numzero; i++)
    A *= ipow(1.0-p[i],data.e[i]); // A = f(1)

  for (int j = 0; j <= data.numcv; j++) {
    for (int k = 1; k < data.d[j]; k++) {
      cdouble R = 0.0;
      for (int i = 0; i <= data.numzero; i++)
	R += cdouble(data.e[i]) / ipow(c[j]-p[i],k);
      y[eq++] = R;
    }
    if (j != data.numcv) {
      cdouble B = 1.0;
      for (int i = 0; i <= data.numzero; i++)
	B *= ipow(c[j]-p[i],data.e[i]);
      y[eq++] = p1distance(B / A, data.v[j]);
    }
  }
}

int solvehd (const hurwitzdata &data, cdouble p[], cdouble c[], int iter, double info[LM_INFO_SZ])
{
  int einf = 0, numcv = 0;
  for (int i = 0; i <= data.numzero; i++)
    if (data.e[i] == 0) // illegal, should all be poles or zeros
      return -1;
    else if (data.e[i] > 0)
      einf += data.e[i], numcv += data.e[i]-1;
  if (einf != data.degree) // degree of numerator should be d
    return -1;

  for (int i = 0; i < data.numzero; i++)
    if (data.e[i] < 0)
      einf += data.e[i], numcv += -data.e[i]-1;
  if (einf <= 0) // order of function at infinity should be at least 1
    return -1;

  for (int i = 0; i <= data.numcv; i++)
    if (data.d[i] < 2) // must all be genuine critical points
      return -1;
    else
      numcv += data.d[i]-1;
  if (numcv+(einf-1) != 2*data.degree-2) // total number of c.p. should be 2d-2
    return -1;

  int m = data.numzero + data.numcv;
  cdouble x0[m];
  for (int i = 0; i < data.numzero; i++)
    x0[i] = p[i];
  for (int i = 0; i < data.numcv; i++)
    x0[i+data.numzero] = c[i];

  iter = dlevmar_dif ((void (*)(double*, double*, int, int, void*)) solvehd_f,
		      (double*) x0, NULL, 2*m, 2*m, iter, NULL, info, NULL, NULL, (void *) &data);

  for (int i = 0; i < data.numzero; i++)
    p[i] = x0[i];
  for (int i = 0; i < data.numcv; i++)
    c[i] = x0[i+data.numzero];

  return iter;
}

int main(void)
{
#if 0 // degree-2 rational map: A*z*(z-p0)/(z-p1)
  const hurwitzdata data = { 2, // degree
			     2, // number of zeros/poles, excluding 0,infty
			     1, // number of critical points, excluding zeros,poles,1
			     { 1, -1, 1 }, // orders of zeros/poles
			     { 2, 2 }, // orders of non-pole/zero critical points
			     { 9.0 } // critical values
  };
  cdouble p[2*MAXDEGREE-2] = { 1.5, 2.0 }; // initial zeros/poles
  cdouble c[2*MAXDEGREE-1] = { 3.1 }; // initial critical points; really 3
#elif 0
 // degree-2 polynomial: A*z*(z-p0)
  const hurwitzdata data = { 2, // degree
			     1, // number of zeros/poles, excluding 0,infty
			     0, // number of critical points, excluding zeros,poles,1
			     { 1, 1 }, // orders of zeros/poles
			     { 2 }, // orders of non-pole/zero critical points
			     { } // critical values
  };
  cdouble p[2*MAXDEGREE-2] = { cdouble(1.0,1.0) }; // initial zeros/poles
  cdouble c[2*MAXDEGREE-2] = { }; // initial critical points
#elif 0 // degree-3 polynomial map: A*z*(z-p0)*(z-p1)
  const hurwitzdata data = { 3, // degree
			     2, // number of zeros/poles, excluding 0,infty
			     1, // number of critical points, excluding zeros,poles,1
			     { 1, 1, 1 }, // orders of zeros/poles
			     { 2, 2 }, // orders of non-pole/zero critical points
			     { -0.411522633 } // critical values; really -100/243
  };
  cdouble p[2*MAXDEGREE-2] = { 4.4, 2.8 }; // initial zeros/poles; really 4,10/2
  cdouble c[2*MAXDEGREE-2] = { 3.1 }; // initial critical points; really 10/3
#elif 0 // degree-3 polynomial map: A*z*(z-p0)*(z-p1)
  const hurwitzdata data = { 3, // degree
			     2, // number of zeros/poles, excluding 0,infty
			     0, // number of critical points, excluding zeros,poles,1
			     { 1, 1, 1 }, // orders of zeros/poles
			     { 3 }, // orders of non-pole/zero critical points
			     { } // critical values
  };
  cdouble p[2*MAXDEGREE-2] = { cdouble(1.5,0.86), cdouble(1.5,-0.87) }; // initial zeros/poles; really (3+-sqrt(-3))/2
  cdouble c[2*MAXDEGREE-2] = { }; // initial critical points
#elif 0 // degree-13 rational map
  const hurwitzdata data = { 13, // degree
			     8, // number of zeros/poles, excluding 0,infty
			     4, // number of critical points, excluding zeros,poles,1
			     { 3, 2, 2, 2, -3, -2, -2, -2, 4 }, // orders of zeros/poles
			     { 3, 2, 2, 2, 4 }, // orders of non-pole/zero critical points
			     { 1.0, 1.0, 1.0, 1.0 } // critical values
  };
  cdouble p[2*MAXDEGREE-2] = {
    // initial zeros
    cdouble(1.1,-1.0),
    cdouble(2.0,-0.2),
    cdouble(0.6,-0.2),
    cdouble(-1.0,-0.8),
    // initial poles
    cdouble(0.5,-0.4),
    cdouble(1.6,-0.5),
    cdouble(0.5,0.0),
    cdouble(-0.6,-0.5)
  };
  cdouble c[2*MAXDEGREE-2] = { // initial critical points
    cdouble(-0.1,-1.0),
    cdouble(0.4,-0.2),
    cdouble(-1.0,-0.2),
    cdouble(2.0,-0.8)
  };
#else
  hurwitzdata data;
  cdouble p[2*MAXDEGREE-2], c[2*MAXDEGREE-2];

  for (;;) {
    char s[1000];
    if (!cin.good())
      strcpy(s,"EOF");
    cin >> s;
    if (!strcmp(s,"END")) {
      break;
    } else if (!strcmp(s,"DEGREE")) {
      cin >> data.degree;
      if (data.degree > MAXDEGREE) {
	cerr << data.degree << " > MAXDEGREE = " << MAXDEGREE << ". Repent.\n";
	return -1;
      }
    } else if (!strcmp(s,"ZEROS/POLES")) {
      cin >> data.numzero;
      for (int i = 0; i <= data.numzero; i++) {
	cin >> data.e[i];
	if (i != data.numzero)
	  cin >> p[i];
      }
    } else if (!strcmp(s,"CRITICAL")) {
      cin >> data.numcv;

      if (data.numcv < 0) { // all critical points are 0,infty, no place for other
	cout << "DEGREE " << data.degree << endl;
	cout << "ZEROS/POLES " << data.numzero << endl;
	for (int i = 0; i <= data.numzero; i++) {
	  cout << data.e[i];
	  if (i != data.numzero) cout << " " << p[i];
	  cout << endl;
	}
	cout << "CRITICAL " << data.numcv << endl;
	cout << "END" << endl;
	return 0;
      }

      for (int i = 0; i <= data.numcv; i++) {
	cin >> data.d[i];
	if (i != data.numcv) {
	  cin >> c[i];
	  cin >> data.v[i];
	}
      }
    } else {
      cerr << "Unknown input " << s << ". Repent." << endl;
      return -1;
    }
  }
#endif

  double info[LM_INFO_SZ];

  int iter = solvehd (data, p, c, maxiter, info);

  if (iter < 0) {
    cerr << "Hurwitz data is invalid. Check again the degrees." << endl;
    return -1;
  }

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

  cout.precision(20);
  cout << "DEGREE " << data.degree << endl;
  cout << "ZEROS/POLES " << data.numzero << endl;
  for (int i = 0; i <= data.numzero; i++) {
    cout << data.e[i];
    if (i != data.numzero) cout << " " << p[i];
    cout << endl;
  }
  cout << "CRITICAL " << data.numcv << endl;
  for (int i = 0; i <= data.numcv; i++) {
    cout << data.d[i];
    if (i != data.numcv) cout << " " << c[i] << " " << data.v[i];
    cout << endl;
  }
  cout << "END" << endl;

  return 0;
}
