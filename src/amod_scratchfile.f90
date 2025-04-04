!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: putrec.f
module mod_scratchfile 
    private
    public :: setrec,putrec,endrec,getrec 
    contains
    subroutine putrec (a,n,nrg)
    implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine writes a matrix to the scratch file on unit 10.
!     -----------------------------------------------------------------
!
    dimension a(n,n)
    irec = nrg+1
    write (unit=10,rec=irec) a
    return
    end

    !     File created at Fri Jun  5 21:58:58 PDT 2020
!     Original source code: setrec.f

    subroutine setrec (n,nnrg)
        implicit double precision (a-h,o-z)
  !
  !     -----------------------------------------------------------------
  !     This subroutine opens a direct access scratch file on unit 10.
  !     -----------------------------------------------------------------
  !
        diskmb = 8.0d-6*n*n*(nnrg+1)
  !     if (diskmb .lt. 1000.d0) then
           open (unit=10,status='scratch',form='unformatted', &
                 access='direct',recl=8*n*n)
           write (6,61) diskmb
  !     else
  !        nnrgmx = (1000.d0/diskmb)*(nnrg+1)-1
  !        write (6,62) nnrgmx
  !        stop
  !     endif
        return
    61  format(/1x,'SETREC:'/1x,70('-')/1x, &
         '*** This run is using ',f8.3,'  Mb of scratch disk space'/1x, &
         70('-'))
    62  format(/1x,'SETREC:'/1x,70('-')/1x, &
         '*** This run needs more than 1 Gb of scratch disk space'/1x, &
         '*** Reduce nnrg to ',i4,' and try again'/1x,70('-'))
        end
    !     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: endrec.f

        subroutine endrec
            implicit double precision (a-h,o-z)
      !
      !     -----------------------------------------------------------------
      !     This subroutine closes and deletes the scratch file on unit 10.
      !     -----------------------------------------------------------------
      !
            close (unit=10,status='delete')
            return
            end
      !     File created at Fri Jun  5 21:58:56 PDT 2020
!     Original source code: getrec.f

            subroutine getrec (a,n,nrg)
                implicit double precision (a-h,o-z)
          !
          !     -----------------------------------------------------------------
          !     This subroutine reads a matrix from the scratch file on unit 10.
          !     -----------------------------------------------------------------
          !
                dimension a(n,n)
                irec = nrg+1
                read (unit=10,rec=irec) a
                return
                end

                
end module mod_scratchfile