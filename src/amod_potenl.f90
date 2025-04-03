!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: potenl.f
module mod_Potenl
      private
      public :: potenl
      contains

      subroutine pes (r12,r23,r31,vev)
            implicit double precision (a-h,o-z)
      !
      !     ----------------------------------------------
      !     User defined PES should be written here in eV.
      !     ----------------------------------------------
      !
            vev = 0.0d0*27.211385d0 ! Your surface goes here.
            return
      end

            
      subroutine potenl (ra,sa,cosa,va,ia)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine returns the potential energy surface va as a
!     function of mass-scaled Jacobi coordinates in arrangement ia.
!     -----------------------------------------------------------------
!
      dimension r(3)
      double precision mass,mtot,mred
      common /masses/ mass(3),mtot,mred
      common /scales/ scale(3),rmlmda
!
      ib = ia+1
      if (ib .gt. 3) ib = 1
      ic = 6-ia-ib
      rap = ra/scale(ia)
      cap = 2.0d0*cosa*rap
      sap = scale(ia)*sa
      sbp = mass(ib)/(mass(ib)+mass(ic))*sap
      scp = mass(ic)/(mass(ib)+mass(ic))*sap
      r(ia) = sap
      r(ib) = sqrt(rap**2-cap*sbp+sbp**2)
      r(ic) = sqrt(rap**2+cap*scp+scp**2)
      call potsub (r,vev)
      va = rmlmda*vev
      return
      end
      !     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: potsub.f

      subroutine potsub (r,vev)
            implicit double precision (a-h,o-z)
      !
      !     -----------------------------------------------------------------
      !     This subroutine chooses which potential to use on the basis
      !     of the atomic masses in common /masses/, and also ensures
      !     that the potential is called with the bond lengths in the
      !     correct order.
      !     -----------------------------------------------------------------
      !
            dimension r(3),s(3),m(3)
            double precision mass,mtot,mred
            common /masses/ mass(3),mtot,mred
      !
            imax = 1
            imin = 1
            do i = 1,3
               m(i) = nint(mass(i))
               if (m(i) .gt. m(imax)) then
                  imax = i
               else if (m(i) .lt. m(imin)) then
                  imin = i
               endif
            enddo
            if (imax .eq. imin) then
               imin = 1
               imid = 2
               imax = 3
            else
               imid = 6-imax-imin
            endif
            s(1) = r(imax)
            s(2) = r(imid)
            s(3) = r(imin)
            mmax = m(imax)
            mmid = m(imid)
            mmin = m(imin)
      !     if (mmin .lt. 1) stop 'potsub 1'
      !     if (mmid .gt. 2) stop 'potsub 2'
      !     if (mmax .le. 2) then
      !        call hh2pot (s,vev)
      !     else if (mmax .eq. 19) then
      !        call fh2pot (s,vev)
      !     else if (mmax .eq. 35) then
      !        call clh2pt (s,vev)
      !     else if (mmax .eq. 37) then
      !        call clh2pt (s,vev)
      !     else
      !        stop 'potsub 3'
      !     endif
            call pes(s(1),s(2),s(3),vev)
            return
            end
      

end module mod_Potenl 