!     File created at Fri Jun  5 21:58:58 PDT 2020
!     Original source code: syminv.f
module mod_syminv 
      private
      public :: syminv
      contains
      subroutine syminv (a,lda,n,ierr)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine uses LAPACK DSYTRF and DSYTRI
!     to invert a real symmetric indefinite matrix.
!     -----------------------------------------------------------------
!
      dimension a(lda,n)

!     dimension ipiv(n),work(32*n)
      allocatable ipiv(:),work(:)

      allocate (ipiv(n),work(32*n))
!
      lwork = 32*n
      call dsytrf ('u',n,a,lda,ipiv,work,lwork,ierr)
      if (ierr .ne. 0) return
      call dsytri ('u',n,a,lda,ipiv,work,ierr)
      do j = 2,n
         do i = 1,j-1
            a(j,i) = a(i,j)
         enddo
      enddo
      return
      end
end module mod_syminv 