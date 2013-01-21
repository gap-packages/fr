/* nag_real_symm_sparse_eigensystem_init (f12fac) Example Program.
 *
 * Copyright 2005 Numerical Algorithms Group.
 *
 * Mark 8, 2005.
 */

#include <nag.h>
#include <nag_stdlib.h>
#include <nag_string.h>
#include <stdio.h>
#include <nagf12.h>
#include <nagf16.h>
static void tv(Integer, double *, double *);
static void av(Integer, double *, double *);

int main(void)
{
  /* Constants */
  Integer licomm=140, imon=0;
  /* Scalars */
  double sigma=0, estnrm;
  Integer exit_status, irevcm, j, lcomm, n, nconv, ncv, nev;
  Integer niter, nshift, nx;
  /* Nag types */
  NagError fail;
  /* Arrays */
  double *comm=0, *eigv=0, *eigest=0;
  double *resid=0, *v=0;
  Integer *icomm=0;
  /* Ponters */
  double *mx=0, *x=0, *y=0;

  exit_status = 0;
  INIT_FAIL(fail);

  Vprintf("nag_real_symm_sparse_eigensystem_init (f12fac) Example Program "
          "Results\n");
  /* Skip heading in data file. */
  Vscanf("%*[^\n] ");

  /* Read values for nx, nev and cnv from data file. */
  Vscanf("%ld%ld%ld%*[^\n] ", &nx, &nev, &ncv);

  /* Allocate memory */
  n = nx * nx;
  lcomm = 3*n + ncv*ncv + 8*ncv + 60;
  if ( !(comm = NAG_ALLOC(lcomm, double)) ||
       !(eigv = NAG_ALLOC(ncv, double)) ||
       !(eigest = NAG_ALLOC(ncv, double)) ||
       !(resid = NAG_ALLOC(n, double)) ||
       !(v = NAG_ALLOC(n * ncv, double)) ||
       !(icomm = NAG_ALLOC(licomm, Integer)) )
    {
      Vprintf("Allocation failure\n");
      exit_status = -1;
      goto END;
    }

  for (j = 0; j < licomm; j++)
   icomm[j] = -9999;
 for (j = 0; j < lcomm; j++)
   comm[j] = -9999.0;

  /* Initialise communication arrays for problem using
     nag_real_symm_sparse_eigensystem_init (f12fac). */
  nag_real_symm_sparse_eigensystem_init(n, nev, ncv, icomm,
                                        licomm, comm, lcomm,
                                        &fail);
  /* Select the required spectrum using
     nag_real_symm_sparse_eigensystem_option (f12fdc). */
  nag_real_symm_sparse_eigensystem_option("smallest magnitude",
                                          icomm, comm, &fail);
  /* Increase the iteration limit if required. */
  nag_real_symm_sparse_eigensystem_option("iteration limit=500",
                                          icomm, comm, &fail);
  irevcm = 0;
 REVCOMLOOP:

 for (j = 0; j < licomm; j++)
   printf("%ld) %ld\n",j,icomm[j]);
 for (j = 0; j < lcomm; j++)
   printf("%ld) %g\n",j,comm[j]);

  /* Repeated calls to reverse communication routine
     nag_real_symm_sparse_eigensystem_iter (f12fbc). */
  nag_real_symm_sparse_eigensystem_iter(&irevcm, resid, v, &x, &y,
                                        &mx, &nshift, comm, icomm,
                                        &fail);
  if (irevcm != 5)
    {
      if (irevcm == -1 || irevcm == 1)
        {
          /* Perform matrix vector multiplication y <--- Op*x */
          av(nx, x, y);
        }
      else if (irevcm == 4 && imon == 1)
        {
          /* If imon=1, get monitoring information using
             nag_real_symm_sparse_eigensystem_monit (f12fec). */
          nag_real_symm_sparse_eigensystem_monit(&niter, &nconv, eigv,
                                                 eigest, icomm, comm);
          /* Compute 2-norm of Ritz estimates using
             nag_dge_norm (f16rac).*/
          nag_dge_norm(Nag_ColMajor, Nag_FrobeniusNorm, nev, 1, eigest,
                       nev, &estnrm, &fail);
          Vprintf("Iteration %3ld, ", niter);
          Vprintf(" No. converged = %3ld,", nconv);
          Vprintf(" norm of estimates = %16.8e\n", estnrm);
        }
      goto REVCOMLOOP;
    }
  if (fail.code == NE_NOERROR)
    {
      /* Post-Process using nag_real_symm_sparse_eigensystem_sol
         (f12fcc) to compute eigenvalues/vectors. */
      nag_real_symm_sparse_eigensystem_sol(&nconv, eigv, v, sigma,
                                           resid, v, comm, icomm,
                                           &fail);
      Vprintf("\n");
      Vprintf("\n The %4ld Ritz values",nconv);
      Vprintf(" of smallest magnitude are:\n\n");
      for (j = 0; j <= nconv-1; ++j)
        {
          Vprintf("%8ld%5s%12.4f\n", j+1, "", eigv[j]);
        }
    }
  else
    {
      Vprintf(" Error from nag_real_symm_sparse_eigensystem_iter (f12fbc)."
              "\n%s\n", fail.message);
      exit_status = 1;
      goto END;
    }
 END:
  if (comm) NAG_FREE(comm);
  if (eigv) NAG_FREE(eigv);
  if (eigest) NAG_FREE(eigest);
  if (resid) NAG_FREE(resid);
  if (v) NAG_FREE(v);
  if (icomm) NAG_FREE(icomm);

  return exit_status;
}

static void av(Integer nx, double *v, double *w)
{
  /* Scalars */
  double nx2;
  Integer j, lo;
  /* Nag types */
  NagError fail;
  /* Function Body */
  INIT_FAIL(fail);
  nx2 = ((double) ((nx + 1) * (nx + 1)));
  tv(nx, v, w);
  nag_daxpby(nx, -nx2, &v[nx], 1, nx2, w, 1, &fail);
  for (j = 1; j <= nx - 2; ++j)
    {
      lo = j * nx;
      tv(nx, &v[lo], &w[lo]);
      nag_daxpby(nx, -nx2, &v[lo-nx], 1, nx2, &w[lo], 1, &fail);
      nag_daxpby(nx, -nx2, &v[lo+nx], 1, 1.0, &w[lo], 1, &fail);
    }
  lo = (nx - 1) * nx;
  tv(nx, &v[lo], &w[lo]);
  nag_daxpby(nx, -nx2, &v[lo-nx], 1, nx2, &w[lo], 1, &fail);
  return;
} /* av */


static void tv(Integer nx, double *x, double *y)
{
  /* Compute the matrix vector multiplication y<---T*x where T is a nx */
  /* by nx tridiagonal matrix with constant diagonals (dd, dl and du). */
  /* Scalars */
  double dd, dl, du;
  Integer j;
  /* Function Body */
  dd = 4.0;
  dl = -1.0;
  du = -1.0;
  y[0] = dd * x[0] + du * x[1];
  for (j = 1; j <= nx - 2; ++j)
    {
      y[j] = dl * x[j-1] + dd * x[j] + du * x[j+1];
    }
  y[nx-1] = dl * x[nx-2] + dd * x[nx-1];
  return;
} /* tv */
