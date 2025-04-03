!     File created at Fri Jun  5 21:58:57 PDT 2020
!     Original source code: rbesjy.f
module Math_Functions
      private
      public :: bessel,rbesjy
      contains
      !     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: bessel.f


      subroutine gensol (a,lda,n,b,ldb,m,ierr)
         implicit double precision (a-h,o-z)
   !
   !     -----------------------------------------------------------------
   !     This subroutine uses LAPACK DGETRF and DGETRS to
   !     solve the linear equations A*X = B.
   !     -----------------------------------------------------------------
   !
         dimension a(lda,n),b(ldb,m)
   
   !     dimension ipiv(n)
         allocatable ipiv(:)
   
         allocate (ipiv(n))
   
   !
         call dgetrf (n,n,a,lda,ipiv,ierr)
         if (ierr .ne. 0) return
         call dgetrs ('N',n,m,a,lda,ipiv,b,ldb,ierr)
         return
         end

      function rgamma(x,odd,even)
         implicit double precision (a-h,o-z)
   !
   !     -----------------------------------------------------------------
   !     Direct fortran translation of Temme's algol routine for computing
   !     rgamma = 1/Gamma(1-x), along with its odd and even parts, for
   !     abs(x) .le. 0.5. [ N.M.Temme, J Comput Phys 19 (1975) 324-337 ]
   !     -----------------------------------------------------------------
   !
         dimension b(12)
         data b / -0.283876542276024d0, -0.076852840844786d0, &
                  +0.001706305071096d0, +0.001271927136655d0, &
                  +0.000076309597586d0, -0.000004971736704d0, &
                  -0.000000865920800d0, -0.000000033126120d0, &
                  +0.000000001745136d0, +0.000000000242310d0, &
                  +0.000000000009161d0, -0.000000000000170d0 /
         save b
   !
         x2 = x*x*8.d0
         alfa = -0.000000000000001d0
         beta = 0.d0
         do i = 12,2,-2
            beta = -(2*alfa+beta)
            alfa = -beta*x2-alfa+b(i)
         enddo
         even = (beta/2.d0+alfa)*x2-alfa+0.921870293650453d0
         alfa = -0.000000000000034d0
         beta = 0.d0
         do i = 11,1,-2
            beta = -(2*alfa+beta)
            alfa = -beta*x2-alfa+b(i)
         enddo
         odd = 2*(alfa+beta)
         rgamma = odd*x+even
         return
         end
      subroutine mbessk (v,x,ck,dk,ek)
         
         implicit double precision (a-h,o-z)
   !
   !     -----------------------------------------------------------------
   !     This subroutine uses Temme's method [ N.M.Temme, J Comput Phys
   !     19 (1975) 324-337 ] to calculate the Modified Bessel function
   !
   !     K(v,x) = ck BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.f9
   !
   !     and its first derivative with respect to x
   !
   !     d/dx K(v,x) = dk BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHP
   !
   !     for a given real order v >= 0 and real argument x > 0.
   !     Note the exponential scaling, which is used to avoid
   !     overflow of K(v,x) for v >> x and underflow for v << x.
   !     -----------------------------------------------------------------
   !
         parameter (eps = 1.d-15)! consistent with rgamma
         parameter (maxit = 1000)
   !
         if (v.lt.0.d0 .or. x.le.0.d0) stop 'mbessk 0'
         pi = acos(-1.d0)
         xmin = 1.d0
   !
   !     begin by calculating K(a,x) and K(a+1,x) for |a| <= 1/2
   !
         na = int(v+0.5d0)
         a = v-na
         if (x .lt. xmin) then
   !
   !        using Temme's series for small x
   !
            b = x/2.d0
            d = -dlog(b)
            e = a*d
            c = a*pi
            if (abs(c) .lt. eps) then
               c = 1.d0
            else
               c = c/sin(c)
            endif
            if (abs(e) .lt. eps) then
               s = 1.d0
            else
               s = sinh(e)/e
            endif
            e = exp(e)
            g = e*rgamma(a,p,q)
            e = (e+1.d0/e)/2.d0
            f = c*(p*e+q*s*d)
            e = a*a
            p = 0.5d0*g*c
            q = 0.5d0/g
            c = 1.d0
            d = b*b
            ak = f
            ak1 = p
            do n = 1,maxit
               f = (f*n+p+q)/(n*n-e)
               c = c*d/n
               p = p/(n-a)
               q = q/(n+a)
               g = c*(p-n*f)
               h = c*f
               ak = ak+h
               ak1 = ak1+g
               if (h/ak+abs(g)/ak1 .lt. eps) go to 1
            enddo
            stop 'mbessk 1'
      1     f = ak
            g = ak1/b
            ex = 0.d0
         else if (x .ge. xmin) then
   !
   !        and Temme's PQ method for large x
   !
            c = 0.25d0-a*a
            g = 1.d0
            f = 0.d0
            e = x*cos(a*pi)/pi/eps
            do n = 1,maxit
               h = (2*(n+x)*g-(n-1+c/n)*f)/(n+1)
               f = g
               g = h
               if (h*n .gt. e) go to 2
            enddo
            stop 'mbessk 2'
      2     p = f/g
            q = p
            b = x+x
            e = b-2.d0
            do m = n,1,-1
               p = (m-1+c/m)/(e+(m+1)*(2.d0-p))
               q = p*(q+1.d0)
            enddo
            f = sqrt(pi/b)/(1.d0+q)
            g = f*(a+x+0.5d0-p)/x
            ex = x
         endif
   !
   !     now recur upwards from K(a,x) to K(v,x),
   !     scaling to avoid overflow along the way
   !
         p = 0.d0
         if (na .gt. 0) then
            y = 2.d0/x
            do n = 1,na
               h = y*(a+n)*g+f
               f = g
               g = h
      3        if (abs(f) .gt. 4.d0) then
                  p = p+1.d0
                  f = 0.0625d0*f
                  g = 0.0625d0*g
                  go to 3
               endif
            enddo
         endif
         ck = f
         dk = (v/x)*f-g
         sk = sqrt(ck*ck+dk*dk)
         ck = ck/sk
         dk = dk/sk
         ek = dlog(sk)+p*dlog(16.d0)-ex
         return
         end

      subroutine rbessk (ell,x,ck,dk,ek)
         implicit double precision (a-h,o-z)
   !
   !     -----------------------------------------------------------------
   !     Modified Riccati-Bessel function of the third kind
   !     and its first derivative with respect to x:
   !
   !     k(ell,x) = ck BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.
   !     d/dx k(ell,x) = dk BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LST
   !     -----------------------------------------------------------------
   !
         if (x.le.0.0d0 .or. ell.lt.-0.5d0) stop 'rbessk 0'
         v = ell+0.5d0
         call mbessk (v,x,ck,dk,ek)
         pi = acos(-1.d0)
         ex = 0.5d0*dlog(pi*x/2.d0)
         dk = dk+ck/(2.d0*x)
         sk = sqrt(ck*ck+dk*dk)
         ck = ck/sk
         dk = dk/sk
         ek = ek+dlog(sk)+ex
         return
         end

      subroutine bessel (x,y,z,cvr,cvi,eint, &
         ilev,jlev,klev,llev,nlev,ered)
         use mod_qvib, only: qvib 
         use mod_pvib, only: pvib 

implicit double precision (a-h,o-z)
double precision llev
!
!     -----------------------------------------------------------------
!     This subroutine forms the asymptotic solution and derivative
!     matrices a,b,c, and d needed for matching the asymptotic Delves
!     coordinate log derivative matrix y onto a reactance matrix k,
!     as in eqs. (116) to (121) of Pack and Parker, and then proceeds
!     to calculate the reactance matrix in the array x.
!     -----------------------------------------------------------------
!
!     common blocks
!
common /arrays/ mro,mvi,nvi,n
common /ranges/ rmin,rmax,smax
common /rotors/ jtot,ipar,jpar,jmax,kmin,kmax
!
!     input arrays
!
dimension x(n,n),y(n,n),z(n,n)
dimension cvr(nvi,n),cvi(nvi,n),eint(n)
dimension ilev(n),jlev(n),klev(n),llev(n),nlev(n)
!
!     local arrays
!
!     dimension wvi(mvi),xvi(mvi),pvi(nvi,2)
allocatable wvi(:),xvi(:),pvi(:,:)

!     dimension fr(n),fvi(n),gvi(n)
allocatable fr(:),fvi(:),gvi(:)

!     dimension fa(n),fb(n),fc(n),fd(n)
allocatable fa(:),fb(:),fc(:),fd(:)

!     dimension a(n,0:nvi-1),b(n,0:nvi-1)
allocatable a(:,:),b(:,:)

!     dimension erow(n),ecol(n),srow(n),scol(n)
allocatable erow(:),ecol(:),srow(:),scol(:)

allocate (wvi(mvi),xvi(mvi),pvi(nvi,2))
allocate (fr(n),fvi(n),gvi(n))
allocate (fa(n),fb(n),fc(n),fd(n))
allocate (a(n,0:nvi-1),b(n,0:nvi-1))
allocate (erow(n),ecol(n),srow(n),scol(n))

!
!     initialisation
!
do nj = 0,nvi-1
do i = 1,n
a(i,nj) = 0.d0
b(i,nj) = 0.d0
enddo
enddo
do j = 1,n
do i = 1,n
x(i,j) = 0.d0
z(i,j) = 0.d0
enddo
enddo
!
!     vibrational quadrature rule
!
pi = acos(-1.d0)
piby2 = 0.5d0*pi
rtrho = sqrt(rmax)
tworho = 2.d0*rmax
smin = 0.d0
tmin = 0.d0
tmax = asin(min(1.d0,smax/rmax))
call qvib (tmin,tmax,mvi,wvi,xvi)
!
!     Delves vibrational quadrature
!
do kvi = 1,mvi
weight = rtrho*wvi(kvi)
ta = xvi(kvi)
cta = cos(ta)
sta = sin(ta)
ra = rmax*cta
sa = rmax*sta
if (sa .lt. smax) then
!
!           Delves vibrational functions
!
call pvib (tmin,ta,tmax,cvr,nvi,n,fr,0)
do i = 1,n
fr(i) = weight*fr(i)
enddo
!
!           Jacobi vibrational functions
!
call pvib (smin,sa,smax,cvi,nvi,n,fvi,0)
call pvib (smin,sa,smax,cvi,nvi,n,gvi,1)
!
!           Jacobi translational (Riccati-Bessel) functions
!
do i = 1,n
ell = llev(i)
psq = ered-eint(i)
if (psq .gt. 0.d0) then
  p = sqrt(psq)
  arg = p*ra
  call rbesjy (ell,arg,cj,dj,ej,cy,dy,ey)
!
!                 exponential scaling of j and y
!                 (for stability near channel thresholds)
!
  if (kvi .eq. 1) then
     ecol(i) = ej
     scol(i) = exp(ej)
     erow(i) = ey
     srow(i) = exp(-ey)
  else
     sj = exp(ej-ecol(i))
     cj = sj*cj
     dj = sj*dj
     sy = exp(ey-erow(i))
     cy = sy*cy
     dy = sy*dy
  endif
  rtp = sqrt(p)
  atr = cj/rtp
  btr = cy/rtp
  ctr = dj*rtp
  dtr = dy*rtp
else
  p = sqrt(-psq)
  arg = p*ra
  call rbessk (ell,arg,ck,dk,ek)
!
!                 exponential scaling of k
!                 (for the same reason)
!
  if (kvi .eq. 1) then
     ecol(i) = 0.d0
     scol(i) = 0.d0
     erow(i) = ek
     srow(i) = 0.d0
  else
     sk = exp(ek-erow(i))
     ck = sk*ck
     dk = sk*dk
  endif
  rtp = sqrt(p)
  atr = 0.d0
  btr = ck/rtp
  ctr = 0.d0
  dtr = dk*rtp
endif
!
!              Pack and Parker eqs. (116) to (121)
!
fa(i) = atr*fvi(i)
fb(i) = btr*fvi(i)
fc(i) = cta*ctr*fvi(i)+sta*atr*gvi(i)
fd(i) = cta*dtr*fvi(i)+sta*btr*gvi(i)
fc(i) = fa(i)/tworho+fc(i)
fd(i) = fb(i)/tworho+fd(i)
enddo
!
!           integral accumulation
!
do j = 1,n
ij = ilev(j)
jj = jlev(j)
kj = klev(j)
nj = nlev(j)
do i = 1,n
  ii = ilev(i)
  ji = jlev(i)
  ki = klev(i)
  if (ii.eq.ij .and. ji.eq.jj .and. ki.eq.kj) then
     a(i,nj) = a(i,nj)+fr(i)*fa(j)
     b(i,nj) = b(i,nj)+fr(i)*fb(j)
     x(i,j)  = x(i,j) -fr(i)*fc(j)
     z(i,j)  = z(i,j) -fr(i)*fd(j)
  endif
enddo
enddo
endif
enddo
!
!     y to k
!
do j = 1,n
ij = ilev(j)
jj = jlev(j)
kj = klev(j)
nj = nlev(j)
do k = 1,n
ik = ilev(k)
jk = jlev(k)
kk = klev(k)
if (ik.eq.ij .and. jk.eq.jj .and. kk.eq.kj) then
do i = 1,n
  x(i,j) = x(i,j)+y(i,k)*a(k,nj)
  z(i,j) = z(i,j)+y(i,k)*b(k,nj)
enddo
endif
enddo
enddo
call gensol (z,n,n,x,n,n,ierr)
if (ierr .ne. 0) stop 'bessel 1'
!
!     elimination of exponential scaling factors
!     from the reactance matrix
!
do j = 1,n
do i = 1,n
x(i,j) = srow(i)*x(i,j)*scol(j)
enddo
enddo
return
end


      subroutine rbesjy (ell,x,cj,dj,ej,cy,dy,ey)
      implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     Riccati-Bessel functions of fractional order
!     and their first derivatives with respect to x:
!
!     j(ell,x) = cj BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.
!     y(ell,x) = cy BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.
!     d/dx j(ell,x) = dj BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LST
!     d/dx y(ell,x) = dy BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LST
!     -----------------------------------------------------------------
!
      if (x.le.0.0d0 .or. ell.lt.-0.5d0) stop 'rbesjy 0'
      v = ell+0.5d0
      call bessjy (v,x,cj,dj,ej,cy,dy,ey)
      pi = acos(-1.d0)
      ex = 0.5d0*dlog(pi*x/2.d0)
      dj = dj+cj/(2.d0*x)
      sj = sqrt(cj*cj+dj*dj)
      cj = cj/sj
      dj = dj/sj
      ej = ej+dlog(sj)+ex
      dy = dy+cy/(2.d0*x)
      sy = sqrt(cy*cy+dy*dy)
      cy = cy/sy
      dy = dy/sy
      ey = ey+dlog(sy)+ex
      return
      end


      !     File created at Fri Jun  5 21:58:55 PDT 2020
!     Original source code: bessjy.f

      subroutine bessjy (v,x,cj,dj,ej,cy,dy,ey)
            implicit double precision (a-h,o-z)
      !
      !     -----------------------------------------------------------------
      !     This subroutine uses a combination of methods (mostly due
      !     to Temme) to calculate the Ordinary Bessel functions
      !
      !     J(v,x) = cj BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.f9
      !     Y(v,x) = cy BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHPOT.f9
      !
      !     and their first derivatives with respect to x
      !
      !     d/dx J(v,x) = dj BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHP
      !     d/dx Y(v,x) = dy BW.3p LSTH.f LSTH.f90 LSTHDATA.f LSTHDATA.f90 LSTHPOT.f LSTHP
      !
      !     for a given real order v >= 0 and real argument x > 0.
      !     Note the exponential scaling, which is used to avoid
      !     overflow of Y(v,x) and underflow of J(v,x) for v >> x.
      !     -----------------------------------------------------------------
      !
            parameter (eps = 1.d-15)! consistent with rgamma
            parameter (maxit = 1000)
      !
            if (v.lt.0.d0 .or. x.le.0.d0) stop 'bessjy 0'
            pi = acos(-1.d0)
            xmin = 3.d0
            xmax = 5.d0-dlog10(eps)
      !
      !     begin by calculating Y(a,x) and Y(a+1,x) for |a| <= 1/2
      !
            na = int(v+0.5d0)
            a = v-na
            if (x .lt. xmin) then
      !
      !        using Temme's series (bessya) for small x
      !        [ N.M.Temme, J Comput Phys 21 (1976) 343-350 ]
      !
               b = x/2.d0
               d = -dlog(b)
               e = a*d
               if (abs(a) .lt. eps) then
                  c = 1.d0/pi
               else
                  c = a/sin(a*pi)
               endif
               if (abs(e) .lt. eps) then
                  s = 1.d0
               else
                  s = sinh(e)/e
               endif
               e = exp(e)
               g = e*rgamma(a,p,q)
               e = (e+1.d0/e)/2.d0
               f = 2*c*(p*e+q*s*d)
               e = a*a
               p = g*c
               q = 1.d0/g/pi
               c = a*pi/2.d0
               if (abs(c) .lt. eps) then
                  r = 1.d0
               else
                  r = sin(c)/c
               endif
               r = pi*c*r*r
               c = 1.d0
               d = -b*b
               ya = f+r*q
               ya1 = p
               do n = 1,maxit
                  f = (f*n+p+q)/(n*n-e)
                  c = c*d/n
                  p = p/(n-a)
                  q = q/(n+a)
                  g = c*(f+r*q)
                  h = c*p-n*g
                  ya = ya+g
                  ya1 = ya1+h
                  del = abs(g)/(1.d0+abs(ya))
                  del1 = abs(h)/(1.d0+abs(ya1))
                  if (del+del1 .lt. eps) go to 1
               enddo
               stop 'bessjy 1'
         1     f = -ya
               g = -ya1/b
            else if (x.ge.xmin .and. x.lt.xmax) then
      !
      !        Temme's PQ method (besspqa) for intermediate x
      !        [ N.M.Temme, J Comput Phys 21 (1976) 343-350 ]
      !
               c = 0.25d0-a*a
               b = x+x
               p = pi
               e = (x*cos(a*pi)/pi/eps)**2
               p = 1.d0
               q = -x
               r = 1.d0+x*x
               s = r
               do n = 2,maxit
                  d = (n-1+c/n)/s
                  p = (2*n-p*d)/(n+1)
                  q = (-b+q*d)/(n+1)
                  s = p*p+q*q
                  r = r*s
                  if (r*n*n .gt. e) go to 2
               enddo
               stop 'bessjy 2'
         2     p = p/s
               f = p
               q = -q/s
               g = q
               do m = n,1,-1
                  r = (m+1)*(2.d0-p)-2.d0
                  s = b+(m+1)*q
                  d = (m-1+c/m)/(r*r+s*s)
                  p = d*r
                  q = d*s
                  e = f+1.d0
                  f = p*e-g*q
                  g = q*e+p*g
               enddo
               f = 1.d0+f
               d = f*f+g*g
               pa = f/d
               qa = -g/d
               d = a+0.5d0-p
               q = q+x
               pa1 = (pa*q-qa*d)/x
               qa1 = (qa*q+pa*d)/x
               b = x-pi*(a+0.5d0)/2.d0
               c = cos(b)
               s = sin(b)
               d = sqrt(2.d0/x/pi)
               f = d*(pa*s+qa*c)
               g = d*(qa1*s-pa1*c)
            else if (x .ge. xmax) then
      !
      !        and Hankel's asymptotic expansions for large x
      !        [ Abramowitz and Stegun, Section 9.2 ]
      !
               p = 0.d0
               q = 0.d0
               do ia = 0,1
                  pa = p
                  qa = q
                  y = 4.d0*(a+ia)**2
                  z = 8.d0*x
                  d = 0.d0
                  w = -1.d0
                  p = 1.d0
                  q = 0.d0
                  tp = 1.d0
                  do k = 1,maxit
                     d = d+z
                     w = w+2.d0
                     tq = +tp*(y-w*w)/d
                     q = q+tq
                     d = d+z
                     w = w+2.d0
                     tp = -tq*(y-w*w)/d
                     p = p+tp
                     if (abs(tp)+abs(tq) .lt. eps) go to 3
                  enddo
                  stop 'bessjy 3'
         3        p = p-0.5d0*tp
                  q = q-0.5d0*tq
               enddo
               pa1 = p
               qa1 = q
               b = x-pi*(a+0.5d0)/2.d0
               c = cos(b)
               s = sin(b)
               d = sqrt(2.d0/x/pi)
               f = d*(pa*s+qa*c)
               g = d*(qa1*s-pa1*c)
            endif
      !
      !     now recur upwards from Y(a,x) to Y(v,x),
      !     scaling to avoid overflow along the way
      !
            p = 0.d0
            if (na .gt. 0) then
               y = 2.d0/x
               do n = 1,na
                  h = y*(a+n)*g-f
                  f = g
                  g = h
         4        if (abs(f) .gt. 4.d0) then
                     p = p+1.d0
                     f = 0.0625d0*f
                     g = 0.0625d0*g
                     go to 4
                  endif
               enddo
            endif
            cy = f
            dy = (v/x)*f-g
            sy = sqrt(cy*cy+dy*dy)
            cy = cy/sy
            dy = dy/sy
            ey = dlog(sy)+p*dlog(16.d0)
      !
      !     finally, calculate J(v,x) and dJ(v,x)/dx
      !
            vv = max(xmin,v)
            if (x .ge. vv) then
      !
      !        using upward recursion in the classically allowed region
      !
               f = d*(pa*c-qa*s)
               g = d*(qa1*c+pa1*s)
               if (na .gt. 0) then
                  y = 2.d0/x
                  do n = 1,na
                     h = y*(a+n)*g-f
                     f = g
                     g = h
                  enddo
               endif
               cj = f
               dj = (v/x)*f-g
               sj = sqrt(cj*cj+dj*dj)
               cj = cj/sj
               dj = dj/sj
               ej = dlog(sj)
            else
      !
      !        and CF1 in the classically forbidden region
      !        [ Numerical Recipes, 2nd Edition, Section 6.7 ]
      !
               ap = 1.d0
               a = v/x
               bp = 0.d0
               b = 1.d0
               f = 0.d0
               g = 0.d0
               y = 2.d0/x
               w = y/pi
               do n = 1,maxit
                  an = y*(v+n)*a-ap
                  ap = a
                  a = an
                  bn = y*(v+n)*b-bp
                  bp = b
                  b = bn
                  if (abs(b) .gt. abs(a)) then
                     ap = ap/b
                     a = a/b
                     bp = bp/b
                     b = 1.d0
                     if (abs(a-f) .lt. eps*abs(f)) then
                        cj = w/(dy-cy*a)
                        dj = a*cj
                        go to 5
                     endif
                     f = a
                  else
                     bp = bp/a
                     b = b/a
                     ap = ap/a
                     a = 1.d0
                     if (abs(b-g) .lt. eps*abs(g)) then
                        dj = w/(dy*b-cy)
                        cj = b*dj
                        go to 5
                     endif
                     g = b
                  endif
               enddo
               stop 'bessjy 4'
         5     sj = sqrt(cj*cj+dj*dj)
               cj = cj/sj
               dj = dj/sj
               ej = dlog(sj)-ey
            endif
            return
            end
      
      END module  