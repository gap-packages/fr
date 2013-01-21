/*     compute eigenvalues of sparse matrix
 *     ----------- input:
 *    N EV CV D
 *    x_1 y_1
 *    ...
 *    x_D y_D
 *    ----------- output:
 *    eigenvalue_1
 *    ...
 *    eigenvalue_EV
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
  Integer licomm=140, lcomm, irevcm, j, n, nconv, ncv, nev,
    nshift, nd, niter, *icomm = 0, *ix = 0, *iy = 0, exit_status;
  double *oldx = 0, *comm = 0, *eigv = 0, *eigest = 0, *mx = 0, *resid = 0,
    *v = 0, *x = 0, *y = 0, sigma;
  NagError fail;

  INIT_FAIL(fail);
  Vscanf("%ld%ld%ld%*[^\n] ", &n, &nev, &ncv, &nd);
  lcomm = 3*n + ncv*ncv + 8*ncv + 60;
  if ( !(oldx = NAG_ALLOC(n, double)) ||
       !(comm = NAG_ALLOC(lcomm, double)) ||
       !(icomm = NAG_ALLOC(licomm, Integer)) ||
       !(eigv = NAG_ALLOC(ncv, double)) ||
       !(eigest = NAG_ALLOC(ncv, double)) ||
       !(mx = NAG_ALLOC(n, double)) ||
       !(resid = NAG_ALLOC(n, double)) ||
       !(v = NAG_ALLOC(n * ncv, double)) ||
       !(x = NAG_ALLOC(n, double)) ||
       !(y = NAG_ALLOC(n, double)) ||
       !(ix = NAG_ALLOC(n, Integer)) ||
       !(iy = NAG_ALLOC(n, Integer)))
    {
      Vprintf("Allocation failure\n");
      exit_status = -1;
      goto END;
    }

  for (j=0; j < nd; j++)
    Vscanf("%d%d", ix+j, iy+j);

  INIT_FAIL(fail);

  f12fac(n, nev, ncv, icomm, licomm, comm, lcomm, &fail);
  f12fdc("largest magnitude", icomm, comm, &fail);
  f12fdc("iteration limit=10000", icomm, comm, &fail);
  irevcm = 0;
 REVCOMLOOP:
  f12fbc(&irevcm, resid, v, &x, &y, &mx, &nshift, comm, icomm, &fail);
  if (irevcm == -1 || irevcm == 1)
    {
      sigma = 0.0;
      for (j=0; j<n; j++) sigma += x[j];
      sigma /= n;
      for (j=0; j<n; j++) { oldx[j] = x[j] - sigma; y[j] = 0.0; }
      for (j=0; j<n; j++) y[ix[j]] = y[ix[j]] + oldx[iy[j]];
      goto REVCOMLOOP;
    }
  if (irevcm == 4)
    {
      f12fec(&niter, &nconv, eigv, eigest, icomm, comm);
      Vprintf("# %d %d\n",niter, nconv);
      goto REVCOMLOOP;
    }
  if (fail.code == NE_NOERROR)
    {
      f12fcc(&nconv, eigv, v, sigma, resid, v, comm, icomm, &fail);
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
  if (oldx) NAG_FREE(oldx);
  if (comm) NAG_FREE(comm);
  if (icomm) NAG_FREE(icomm);
  if (eigv) NAG_FREE(eigv);
  if (eigest) NAG_FREE(eigest);
  if (mx) NAG_FREE(mx);
  if (resid) NAG_FREE(resid);
  if (v) NAG_FREE(v);
  if (x) NAG_FREE(x);
  if (y) NAG_FREE(y);
  if (ix) NAG_FREE(ix);
  if (iy) NAG_FREE(iy);

  return exit_status;
}
