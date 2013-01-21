      program f12fafe
!     compute eigenvalues of sparse matrix
!     ----------- input:
!     N EV CV D
!     x_1 y_1
!     ...
!     x_D y_D
!     ----------- output:
!     eigenvalue_1
!     ...
!     eigenvalue_EV
      implicit none

      integer :: licomm=140, lcomm, ldv
      integer, allocatable :: icomm(:), ix(:), iy(:)
      integer ifail, irevcm, j, n, nconv, ncv, nev, nshift, nd, niter
      double precision, allocatable :: oldx(:), comm(:), d(:,:),
     + mx(:), resid(:), v(:,:), x(:)
      double precision sigma
      external f12faf, f12fbf, f12fcf, f12fdf, f12fef

      read *, n, nev, ncv, nd

      lcomm = 3*n+ncv*ncv+8*ncv+60
      ldv = n

      allocate (oldx(1:n), comm(1:lcomm), icomm(1:licomm), d(1:ncv,1:2),
     + mx(1:n), resid(1:n), v(1:ldv,1:ncv), x(1:n),ix(1:nd), iy(1:nd))
      do j=1,nd
         read *, ix(j), iy(j)
      end do
      ifail = 0
      call f12faf(n,nev,ncv,icomm,licomm,comm,lcomm,ifail)
      call f12fdf('largest magnitude',icomm,comm,ifail)
      call f12fdf('iteration limit=10000',icomm,comm,ifail)
      irevcm = 0
      ifail = -1
 10   call f12fbf(irevcm,resid,v,ldv,x,mx,nshift,comm,icomm,ifail)
      if (ifail.gt.0) then
         print *, 'Error ',ifail,' occurred'
         stop
      end if
      if ((irevcm.eq.-1).or.(irevcm.eq.1)) then
!         print *,'# input vector: ',x(1:n)
         oldx(1:n) = x(1:n) - sum(x)/n
         x(1:n) = 0.d0
         do j=1,nd
            x(ix(j)) = x(ix(j)) + oldx(iy(j))
         end do
!         print *,'# output vector: ',x(1:n)
         goto 10
      end if
      if (irevcm.eq.4) then
         call f12fef(niter,nconv,d,d(1,2),icomm,comm)
         print *, '# ',niter,nconv
         goto 10
      end if
      if (irevcm.eq.5) then
         call f12fcf(nconv,d,v,ldv,sigma,resid,v,ldv,comm,icomm,ifail)
         do j=1,nconv
            print *, '# eigenvector ',j,': ',v(1:n,j)
         end do
         print *, '\\\\\\Number of data points: ',nconv
         do j=1,nconv
            print *, d(j,1)
         end do
         stop
      end if
      print *, 'Unknown case irevcm = ',irevcm
      stop
      end
