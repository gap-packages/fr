      implicit none
      integer maxn, lwork, lrwork
      parameter (maxn=1030, lwork=maxn*(maxn+3), lrwork=1)
      integer i, n, m, k, noits, novecs, liwork, iwork(6*maxn), ifail
      real x(maxn,maxn), d(maxn), work(lwork), rwork(lrwork)
      real tol
      real dot
      external dot, image, monit, f02fjf

      read *, m
      read *, liwork
      do i=1,liwork
         read *, iwork(2*i-1), iwork(2*i)
      end do
      liwork = liwork*2
      n = 0
      do i=1,liwork
         if (iwork(i) .gt. n) n = iwork(i)
      end do
      k = m+4
      if (k .gt. n) k = n
      noits = 200
      tol = 1.e-2
      novecs = 0
      ifail = 1

      print *, 'Running f02fjf with n=',n,' m=',m,' k=',k,' tol=',tol

      call f02fjf(n,m,k,noits,tol,dot,image,monit,novecs,x,maxn,d,work,
     +     lwork,rwork,lrwork,iwork,liwork,ifail)
      if (ifail .ne. 0) then
         print *, 'f02fjf failed with error ', ifail
*         stop
      end if
      do i=1,m
         print *, d(i)
      end do
      stop
      end

      subroutine test(n,iwork,liwork)
      implicit none
      integer n, liwork, lrwork, iwork(liwork)
      real w(10), z(10), rwork(1)
      integer i, j
      do i=1,n
         do j=1,n
            z(j) = 0.e0
         end do
         z(i) = 1.e0
         lrwork = 1
         call image(j,n,z,w,rwork,lrwork,iwork,liwork)
         do j=1,n
            print *,w(j)
         end do
         print *, '---'
      end do
      return
      end

      real function dot(iflag,n,z,w,rwork,lrwork,iwork,liwork)
      implicit none
      integer i, iflag, liwork, lrwork, n, iwork(liwork)
      real rwork(lrwork), w(n), z(n)
      real s
      s = 0.e0
      do i=2,liwork,2
         s = s + z(iwork(i-1))*w(iwork(i))
      end do
      dot = s
      return
      end

      subroutine image(iflag,n,z,w,rwork,lrwork,iwork,liwork)
      implicit none
      integer i, iflag, liwork, lrwork, n, iwork(liwork)
      real rwork(lrwork), w(n), z(n)
      do i=1,n
         w(i) = 0.e0
      end do
      do i=2,liwork,2
         w(iwork(i)) = w(iwork(i)) + z(iwork(i-1))
      end do
      return
      end

      subroutine monit(istate,nextit,nevals,nevecs,k,f,d)
      implicit none
      integer istate, nextit, nevals, nevecs, k
      real f(k), d(k)

      print *, 'Monitor: ', istate, nextit, nevals, nevecs
      return
      end
