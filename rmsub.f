!    part of the frt distribution
c   mathematical helping functions for the frt model
c     no header files included in the subroutines
c
      function eint(xy)
c   exponential integral
c
      implicit none
      double precision eps, eint, xy
c
      data eps/.1d-6/
c
      if (abs(xy) .lt. eps) then
         eint = 1.d0 - xy*(.5d0 - xy/6.d0*(1.d0 - xy*.25d0))
      else
         eint = (1.d0 - exp(-xy))/xy
      endif
      return
      end
*
      subroutine gmfres (calp2, rn, pkhair, gmf)
c  Fresnel' reflection                    A. Kuusk 02.01.1991
c  input parameters are calp2 = cos(th_incident),  rn = refract_ind.,
c  pkhair = leaf hair index
c
c
      implicit none
      double precision calp2, rn, pkhair, gmf
      double precision pi, pi12, x2, ag, bg, xy, cg, sa2, yg, yy
c
      pi = acos(-1.d0)
      pi12 = pi/2.d0
      x2  = calp2*calp2
      ag  = x2*2.d0 - 1.d0 + rn*rn
      bg  = 1.d0 + (ag - 2.d0)*x2
      xy  = ag - x2
      cg  = 2.d0*calp2*sqrt(xy)
      sa2 = 1.d0 - x2
      yg  = (bg + sa2*cg)*(ag + cg)
      yg  = (ag - cg)*bg/yg
      yy  = sqrt(sa2)/pi12/calp2*pkhair
      gmf = exp(-yy)*yg
c
      return
      end
*
      subroutine iterats(clmpx, fgx)
c    input: clmpx -- clumping index c_B
c    output: fgx -- Fischer's grouping index
C*************************************************************)
C*      ITERATS.FOR  on Newton-Raphson meetodil              *)
C*        funktsiooni  Ci-lnGi/(1-Gi)=0 va"a"rtuste Gi       *)
C*               interpoolimiseks  (22.04.97)  Anne Jo~eveer *)
C*************************************************************)
c
      implicit double precision (a-h, o-z)
      parameter (nfgi = 15, maxit=100)
      save clmpj
c
      dimension fgi(1:nfgi),  clmpj(1:nfgi)
c
      data fgi/0.001d0, 0.05d0, 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.6d0,
     &        0.8d0, 1.d0, 1.2d0, 1.5d0, 2.d0, 2.5d0, 3.5d0, 6.d0/
c
c
         do j = 1, nfgi
            clmpj(j) = 1.d0
            if (fgi(j) .ne. 1.d0)
     &         clmpj(j) = -log(fgi(j))/(1.d0 - fgi(j))
         enddo
c
      fgx    = 1.d0
      if (abs(1.d0 - clmpx) .lt. 0.01d0) goto 18
c
      g1     = 0.001d0
      do i = 1, 14
         if (clmpx .le. clmpj(i) .and. clmpx .gt. clmpj(i+1)) then
            g1   = fgi(i)
            g2   = fgi(i+1)
            xacc = 0.1d-5
            call rtsafe(clmpx, g1, g2, xacc, maxit, fgx)
            goto 18
         endif
      enddo
      if (clmpx .lt. clmpj(nfgi)) fgx = 6.d0
c
18    continue
      end
*
      subroutine rtsafe(xxin, xx1, xx2, xacc, maxit, yyout)
c
c  Newton-Raphson method. Numerical Recipes, 9.4
*        Using a combination of Newton-Raphson and bisection,        *)
*        find the root of a function bracketed between xx1 and xx2.  *)
*        The root, returned as yyout                                 *)
*        will be refined until its accuracy is known within +-xacc.  *)
*        funcd is a user-supplied subroutine which returns both the  *)
*        function value and the first derivative of the function.    *)
      implicit none
      double precision xxin, xx1, xx2, xacc, yyout
      integer maxit
      
      double precision fl, df, fh, xl, xh, dx, dxold,
     & f, temp
      integer j
      external funcd
c
      call funcd(xxin, xx1, fl, df)
      call funcd(xxin, xx2, fh, df)
c
      if (fl .eq. 0.d0) then
         yyout = xx1
         return
      else if (fh .eq. 0.d0) then
         yyout = xx2
         return
      else if (fl .lt. 0.d0) then
*                             ! Orient the search so that f(xl) < 0. *)
         xl = xx1
         xh = xx2
      else
         xh = xx1
         xl = xx2
      endif
      yyout  = .5d0*(xx1 + xx2)
      dxold  = abs(xx2 - xx1)
      dx     = dxold
*                             ! Initialize the guess for root,  *)
*                             ! the "stepsize before last,"     *)
*                             !  and the last step.             *)
c
      call funcd(xxin, yyout, f, df)
c
      do 11 j = 1, maxit
*                             !  Loop over allowed iterations.   *)
        if (((yyout - xh)*df - f)*((yyout - xl)*df - f) .ge. 0.d0
     &  .or. abs(2.d0*f) .gt. abs(dxold*df)) then
            dxold = dx
*                             ! Bisect if Newton out of range,   *)
          dx = 0.5d0*(xh - xl)
*                             ! or not decreasing fast enough.   *)
          yyout = xl + dx
          if (xl .eq. yyout) return
*                             ! Change in root is negligible. *)
        else
*                             ! Newton step acceptable. Take it. *)
           dxold = dx
           dx    = f/df
           temp  = yyout
           yyout = yyout - dx
           if (temp .eq. yyout) return
        endif
        if (abs(dx).lt.xacc) return
c
        call funcd(xxin, yyout, f, df)
c
*                             !Convergence criterion.       *)
*                             !The one new function evalution per iteration.*)
*                             !Maintain the bracket on the root *)
        if (f .lt. 0.d0) then
           xl = yyout
        else
           xh = yyout
        endif
11    continue
      return
      end
*
      subroutine funcd(clmpj, fgi, fcg, dfcg)
c     function for iteratively finding the grouping index.
c     called by rtsafe
*     Ci+lnGi/(1-Gi)=0......f(x)=0,  Gi...x  *)
      implicit double precision (a-h, o-z)
c
      xln  = log(fgi)
      vahe = 1.d0 - fgi
      if (abs(vahe) .lt. 0.0001d0) then
        return
      else
        fcg  = clmpj + xln/vahe
        dfcg = (1.d0/fgi - 1.d0 - xln)/vahe/vahe
      endif
      return
      end
*
      subroutine gauleg(x1, x2, x, w, n)
c   Gauss-Legendre n-point quadrature. Press et al., 1992, Chpt. 4.5
c   frt.f: call gauleg(x1, x2, zctst, acztst, ntstz)
!f2py intent(in) x1, x2, n
!f2py intent(out) x, w
! ?? x=points, w=weights (guess)
c
      implicit none
      integer n
      double precision x1, x2, x(n), w(n)
      
      double precision eps, pi, xm, xl, z, z1
      double precision p1, p2, p3, pp
      integer m, i, j
      parameter (eps = .1d-14)
c
      pi = acos(-1.d0)
      m  = (n + 1)/2
      xm = 0.5d0*(x2 + x1)
      xl = 0.5d0*(x2 - x1)
      do i = 1, m
        z = cos(pi*(dble(i) - .25d0)/(dble(n) + .5d0))
10      continue
          p1 = 1.d0
          p2 = 0.d0
          do j = 1, n
            p3 = p2
            p2 = p1
            p1 = ((2.d0*dble(j) - 1.d0)*z*p2 -
     &           (dble(j) - 1.d0)*p3)/dble(j)
          enddo
          pp = dble(n)*(z*p1 - p2)/(z*z - 1.d0)
          z1 = z
          z  = z1 - p1/pp
        if(abs(z - z1) .gt. eps) goto 10
        if (abs(z) .lt. eps) z = 0.d0
        x(i)     = xm - xl*z
        x(n+1-i) = xm + xl*z
        w(i)     = 2.d0*xl/((1.d0 - z*z)*pp*pp)
        w(n+1-i) = w(i)
      enddo
      return
      end
      
      subroutine cubell9(ntst, xtst, ytst, ztst, atst)
c
c   knots and weights of the cubature in a sphere
c                                      A. Kuusk  08.01.2003
c   Mysovskih,I.P. Interpolyatsionnye kubaturnye formuly, 
c   Nauka, Moskva 1981, Chpt. 5, Algorithm 16.31.
c
      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c
*     parameter (ncub = 48, nface = 20)
      parameter (nface = 20)
      dimension xtst(1:ncub), ytst(1:ncub), ztst(1:ncub), atst(1:ncub)
      dimension xi(1:ncub), yi(1:ncub), zi(1:ncub)
      dimension i1(1:nface), i2(1:nface), i3(1:nface)
      data 
     & i1/2,  2,  2,  2,  2,  3,  3,  3,  3,  3, 
     &    4,  4,  4,  5,  5,  6,  6,  7,  7,  8/, 
     & i2/4,  4,  5,  6,  7,  9,  9, 10, 11, 12, 
     &    5,  8,  9,  6,  9,  7, 10,  8, 11, 12/,
     & i3/5,  8,  6,  7,  8, 10, 13, 11, 12, 13,
     &    9, 13, 13, 10, 10, 11, 11, 12, 12, 13/
c
*     
c      pi = acos(-1.d0)
      ntst = 45
      vol  = pi*4.d0/3.d0
c
      atst(1) = vol*2096.d0/42525.d0
      do itst = 2, 13
         atst(itst) = vol*(491691.d0 + 54101.d0*sqrt(31.d0))/21.0924d6
      enddo
      do itst = 14, 25
         atst(itst) = vol*(491691.d0 - 54101.d0*sqrt(31.d0))/21.0924d6
      enddo
      do itst = 26, 45
         atst(itst) = vol*1331.d0/68.04d3
      enddo
c
      alph = sqrt(81.d0 - 6.d0*sqrt(31.d0))/11.d0
      beta = sqrt(81.d0 + 6.d0*sqrt(31.d0))/11.d0
      gamm = 3.d0/sqrt(11.d0)
c
c                                   *****    ceneter of sphere
      xtst(1) = 0.d0
      ytst(1) = 0.d0
      ztst(1) = 0.d0
c                                   *****    verteces of icosahedron (12)
      xi(2) = 0.d0
      xi(3) = 0.d0
      yi(2) = 0.d0
      yi(3) = 0.d0
      zi(2) = 1.d0
      zi(3) = -1.d0
      do i = 4, 8
         xi(i) = cos(dble(i - 4)*2.d0*pi/5.d0)*2.d0/sqrt(5.d0)
         xi(i+5) = cos(dble(2*(i - 4) + 1)*
     &             2.d0*pi/5.d0)*2.d0/sqrt(5.d0)
         yi(i) = sin(dble(i - 4)*2.d0*pi/5.d0)*2.d0/sqrt(5.d0)
         yi(i+5) = sin(dble(2*(i - 4) + 1)*
     &             2.d0*pi/5.d0)*2.d0/sqrt(5.d0)
         zi(i) = 1.d0/sqrt(5.d0)
         zi(i+5) = -1.d0/sqrt(5.d0)
      enddo
      do i = 2, 13
         xtst(i) = xi(i)*alph
         ytst(i) = yi(i)*alph
         ztst(i) = zi(i)*alph
         xtst(i+12) = xi(i)*beta
         ytst(i+12) = yi(i)*beta
         ztst(i+12) = zi(i)*beta
      enddo
c
c                                  ***** projections of facets' centers
c
      i  = 1
      xz = xi(i1(i)) + xi(i2(i)) + xi(i3(i))
      yz = yi(i1(i)) + yi(i2(i)) + yi(i3(i))
      zz = zi(i1(i)) + zi(i2(i)) + zi(i3(i))
      rz = sqrt(xz**2 + yz**2 + zz**2)
      xtst(26) = xz/rz*gamm
      ytst(26) = yz/rz*gamm
      ztst(26) = zz/rz*gamm
      do i = 2, nface
         xtst(i+25) = (xi(i1(i)) + xi(i2(i)) + xi(i3(i)))/rz*gamm
         ytst(i+25) = (yi(i1(i)) + yi(i2(i)) + yi(i3(i)))/rz*gamm
         ztst(i+25) = (zi(i1(i)) + zi(i2(i)) + zi(i3(i)))/rz*gamm
      enddo
c
*     do j = 1, ntst
*        write (*, '(i2, 4g12.5)')
*    &   j, xtst(j), ytst(j), ztst(j), atst(j)
*     enddo
c
      return
      end

      subroutine cubcirc(ntst, xtst, ytst, atst)
c
c   knots and weights of the cubature in a circle
c                                      A. Kuusk  27.08. 2002
c   Vysotskikh, I.P., Nauka, M., 1981, Ch. 5, P-t. 16.33.
c
      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c
*     parameter (ncub = 48)
      dimension xtst(1:ncub), ytst(1:ncub), atst(1:ncub)
      dimension ri(1:8)
      data ri/1.d0, 2.d0, 4.d0, 5.d0, 7.d0, 8.d0, 10.d0, 11.d0/
c
*     pi = acos(-1.d0)
      ntst = 28
      r1 = sqrt((10.d0 - sqrt(10.d0))/15.d0)
      r2 = sqrt((10.d0 + sqrt(10.d0))/15.d0)
      r3 = sqrt((31.d0 - sqrt(601.d0))/60.d0)
      r4 = sqrt(3.d0/5.d0)
      r5 = sqrt((31.d0 + sqrt(601.d0))/60.d0)
      c1 = (340.d0 + 25.d0*sqrt(10.d0))*pi/10368.d0
      c2 = (340.d0 - 25.d0*sqrt(10.d0))*pi/10368.d0
      c3 = (857.d0 + 12707.d0/sqrt(601.d0))*pi/20736.d0
      c4 = 125.d0*pi/3456.d0
      c5 = (857.d0 - 12707.d0/sqrt(601.d0))*pi/20736.d0
c
      do k = 1, 8
         xtst(k)   = r1*cos(pi*ri(k)/6.d0)
         ytst(k)   = r1*sin(pi*ri(k)/6.d0)
         atst(k)   = c1
         xtst(8+k) = r2*cos(pi*ri(k)/6.d0)
         ytst(8+k) = r2*sin(pi*ri(k)/6.d0)
         atst(8+k) = c2
      enddo
      do k = 1, 4
         xtst(16+k) = r3*cos(pi*dble(k)/2.d0)
         ytst(16+k) = r3*sin(pi*dble(k)/2.d0)
         atst(16+k) = c3
         xtst(20+k) = r4*cos(pi*dble(k)/2.d0)
         ytst(20+k) = r4*sin(pi*dble(k)/2.d0)
         atst(20+k) = c4
         xtst(24+k) = r5*cos(pi*dble(k)/2.d0)
         ytst(24+k) = r5*sin(pi*dble(k)/2.d0)
         atst(24+k) = c5
      enddo
c
      return
      end
      
      subroutine rcone
     & (xi, yj, zk, sthrl, cthrl, sphrl, cphrl, 
     &  aell, vhelko, vhcyl, rlout)
c
c   kaugus p-tist (xi, yj, zk) võra pinnani suunas thet, phi
c   A. Kuusk, 7.05.1998
c   vhelko - koonilise osa ko~rgus,  aell - võra raadius,
c   phi on nurk x-teljest, sphrl = sin(phi), cphrl = cos(phi)
c   sthrl = sin(thet), cthrl = cos(thet)
c
      implicit double precision (a-h, o-z)
      data eps/.1d-6/, eps2/.001d0/
c
      z3 = vhelko + eps2
      z0 = zk - eps2
c
      if (cthrl .ge. 0.d0) then
c                                             ülemine poolsfäär
      if (sthrl .lt. eps) then
c                                                kiir on vertikaalne
         rlout = vhelko - vhelko/aell*sqrt(xi**2 + yj**2) - zk
      else
         ae2 = aell**2
         hh2 = vhelko**2
         if (zk .ge. 0.d0) then
c                                                (x,y,z) on koonuses
            aaa = ae2/hh2*(cthrl**2) - sthrl**2
            bbb = 2.d0*ae2*cthrl/vhelko*(zk/vhelko - 1.d0) 
     &            - 2.d0*sthrl*(xi*cphrl + yj*sphrl)
            ccc = (aell*(1.d0 - zk/vhelko))**2 - (xi**2 + yj**2)
            det = bbb**2 - 4.d0*aaa*ccc
            if (abs(det) .lt. eps) then
               rlout = -bbb*.5d0/aaa
            else
               rlout = (sqrt(det) - bbb)*.5d0/aaa
               z2    = rlout*cthrl + zk
               if (z2 .lt. z0 .or. z2 .gt. z3) 
     &            rlout = (-sqrt(det) - bbb)*.5d0/aaa
            endif
c
         else
c                                                (x,y,z) on silindris
            rl1 = 0.d0
            aaa = sthrl**2
            bbb = 2.d0*sthrl*(xi*cphrl + yj*sphrl)
            ccc = xi**2 + yj**2 - ae2
            det = bbb**2 - 4.d0*aaa*ccc
            if (abs(det) .lt. eps) then
               rl1 = -bbb*.5d0/aaa
            else
               rl1 = (sqrt(det) - bbb)*.5d0/aaa
               z2  = rl1*cthrl + zk
               if (z2 .lt. z0) 
     &            rl1 = (-sqrt(det) - bbb)*.5d0/aaa
            endif
c
            z2  = rl1*cthrl + zk
            if (z2 .le. 0.d0) then
c                                   kiir väljub läbi silindri külgpinna
               rlout = rl1
            else
c                                   kiir väljub läbi koonuse
               rl2 = 0.d0
               rl1 = abs(zk)/cthrl
               x3  = xi + rl1*sthrl*cphrl
               y3  = yj + rl1*sthrl*sphrl
               aaa = ae2/hh2*(cthrl**2) - sthrl**2
               bbb = -2.d0*ae2*cthrl/vhelko 
     &               - 2.d0*sthrl*(x3*cphrl + y3*sphrl)
               ccc = ae2 - (x3**2 + y3**2)
               det = bbb**2 - 4.d0*aaa*ccc
               if (abs(det) .lt. eps) then
                  rl2 = -bbb*.5d0/aaa
               else
                  rl2 = (sqrt(det) - bbb)*.5d0/aaa
                  z2  = rl2*cthrl
                  if (z2 .lt. 0.d0 .or. z2 .gt. z3) 
     &            rl2 = (-sqrt(det) - bbb)*.5d0/aaa
               endif
               rlout = rl1 + rl2
            endif
         endif
      endif
      else
c                                     alumine poolsfäär
      if (sthrl .lt. eps) then
c                                                kiir on vertikaalne
         rlout = vhcyl + zk
      else
         ae2 = aell**2 + eps
         if (zk .ge. 0.d0) then
c                                                (x,y,z) on koonuses
            rl1 = -zk/cthrl
            x2  = xi + rl1*sthrl*cphrl
            y2  = yj + rl1*sthrl*sphrl
            if ((x2**2 + y2**2) .ge. ae2) then
c                                      kiir väljub läbi koonuse pinna
               hh2 = vhelko**2
               aaa = aell**2/hh2*(cthrl**2) - sthrl**2
               bbb = 2.d0*aell**2*cthrl/vhelko*(zk/vhelko - 1.d0) 
     &               - 2.d0*sthrl*(xi*cphrl + yj*sphrl)
               ccc = (aell*(1.d0 - zk/vhelko))**2 - (xi**2 + yj**2)
               det = bbb**2 - 4.d0*aaa*ccc
               if (abs(det) .lt. eps) then
                  rlout = -bbb*.5d0/aaa
               else
                  rlout = (sqrt(det) - bbb)*.5d0/aaa
*                 z2    = rlout*cthrl + zk
*                 if (z2 .gt. zk .or. z2 .lt. 0.d0) 
*                 if (z2 .gt. zk) 
                  if (rlout .lt. 0.d0) 
     &               rlout = (-sqrt(det) - bbb)*.5d0/aaa
               endif
            else
c                                     kiir siseneb silindrisse
               rl1 = -(vhcyl + zk)/cthrl
               x2  = xi + rl1*sthrl*cphrl
               y2  = yj + rl1*sthrl*sphrl
               if ((x2**2 + y2**2) .le. ae2) then
c                                    kiir väljub läbi silindri põhja
                  rlout = rl1
               else
c                                    kiir väljub läbi silindri seina
                  aaa = sthrl**2
                  bbb = 2.d0*sthrl*(xi*cphrl + yj*sphrl)
                  ccc = xi**2 + yj**2 - aell**2
                  det = bbb**2 - 4.d0*aaa*ccc
                  if (abs(det) .lt. eps) then
                     rl1 = -bbb*.5d0/aaa
                  else
                     rl1 = (sqrt(det) - bbb)*.5d0/aaa
                  endif
                  z2  = rl1*cthrl + zk
                  if (z2 .gt. 0.d0)
     &               rl1 = (-sqrt(det) - bbb)*.5d0/aaa
                  rlout = rl1
               endif
            endif
c
         else
c                                             (x,y,z) on silindris
            rl1 = -(vhcyl + zk)/cthrl
            x2  = xi + rl1*sthrl*cphrl
            y2  = yj + rl1*sthrl*sphrl
            if ((x2**2 + y2**2) .le. ae2) then
c                                    kiir väljub läbi silindri põhja
               rlout = rl1
            else
c                                    kiir väljub läbi silindri seina
               aaa = sthrl**2
               bbb = 2.d0*sthrl*(xi*cphrl + yj*sphrl)
               ccc = xi**2 + yj**2 - aell**2
               det = bbb**2 - 4.d0*aaa*ccc
               if (abs(det) .lt. eps) then
                  rl1 = -bbb*.5d0/aaa
               else
                  rl1 = (sqrt(det) - bbb)*.5d0/aaa
               endif
               z2  = rl1*cthrl + zk
               if (z2 .gt. zk + eps)
     &            rl1 = (-sqrt(det) - bbb)*.5d0/aaa
               rlout = rl1
            endif
c
         endif
      endif
      endif
c
*     if (rlout .lt. -.01d0) then
*        print *, 'rlout = ', rlout
*        pause
*     endif
      if (rlout .lt. 0.d0) rlout = 0.d0
      return
      end

      subroutine rlips
     & (xi, yj, zk, sthrl, cthrl, sphrl, cphrl, aell, cell, rlout)
c
c    rlips   A. Kuusk   29.12.1985
c    kaugus p - tist xi, yj, zk ellipsoidi pinnani suunas thet, phi;
c    sthrl = sin(thet), cthrl = cos(thet), 0 <= thet <= pi/2
c    sphrl = sin(phi),  cphrl = cos(phi),  0 <= phi  <= pi
c    phi on nurk x-telje ja kiirt (thet, phi) la"biva vertikaaltasandi vahel.
c    aell, cell - ellipsoidi poolteljed.
c    Koordinaatide algus on ellipsoidi tsentris.
c    
      implicit double precision (a-h, o-z)
c
      data eps/.1d-4/
c
      a2    = aell*aell
      c2    = cell*cell
      aaa   = c2*(sthrl**2) + a2*(cthrl**2)
      bbb   = 2.d0*(c2*sthrl*(xi*cphrl + yj*sphrl) + a2*cthrl*zk)
      ccc   = c2*(xi**2 + yj**2) + a2*(zk**2 - c2)
      det   = (bbb**2 - 4.d0*aaa*ccc)
      if (abs(det) .lt. eps) then
         rlout = -bbb*.5d0/aaa
      else
         rlout = (sqrt(bbb**2 - 4.d0*aaa*ccc) - bbb)*.5d0/aaa
      endif
c
      if (rlout .lt. 0d0) rlout = 0.d0
c
      return
      end
*
      subroutine scone(hcone, hcz, rcone, tgthx, scz)
c
c  (Tüvi)koonuse projektsioon                 A. Kuusk   14.11.2000
c  hcone - koonuse kõrgus, hcz - tüvikoonuse kõrgus, rcone - koonuse
c  aluse raadius, tgthx - vaatenurga tangens
c
      implicit double precision (a-h, o-z)
      data eps/.1d-4/

      pi = acos(-1.d0)
      rhc   = rcone/hcone/max(eps, tgthx)
c                                              thx > thc
      if (rhc .lt. 1.d0) then
         hc4  = abs(hcone - hcz)
         beta = acos(rcone/(hcone*tgthx))
c                                        tüvikoonus
         if (hc4 .gt. eps) then
            rcz = rcone*hc4/hcone
            sc3 = rcz**2*beta
            sc4 = sqrt((hc4*tgthx)**2 - rcz**2)*rcz
         else
            sc3 = 0.d0
            sc4 = 0.d0
         endif
c
         sc1 = sqrt((hcone*tgthx)**2 - rcone**2)*rcone
         sc2 = rcone**2*(pi - beta)
         scz = sc2 + sc1 - sc4 + sc3
      else
         scz  = pi*rcone**2 
      endif
      return
      end
******************************************************************** scone
*
      subroutine stem(z1, z2, dbh, htree, styvi)
c
c  Tüve pikilõike pindala                A. Kuusk   20.04.2005
c  dbh - rinnasdiameeter, htree - puu kõrgus
c  Pindala arvutatakse nivoost z1 nivooni z2 ristlõike 
c  diameetri integreerimisega. Tüve moodustaja valemid on
c  pärit lätlaste 1988. a. normatiivdokumentidest.
c
      implicit double precision (a-h, o-z)
      include 'frtpar.h'
      dimension ati(1:7)
c  need parameetrid on männi jaoks
      data ati/118.981d0, -277.578d0, 1140.525d0, 
     & -3037.487d0, 4419.682d0, -3361.78d0, 997.657d0/
      data  ht0/26.d0/, dt0/30.d0/, pty/0.007d0/, qty/-0.007d0/
c
      xz0 = 1.3d0/htree
      xz1 = z1/htree
      xz2 = z2/htree
      eet = pty*(htree - ht0) + qty*(dbh - dt0)
c
      sum1 = 0.d0
      do j = 1, 7
         sum1 = sum1 + ati(j)*(xz0**(j - 1))
      enddo
      f13 = sum1*(1.d0 + eet*(xz0*xz0 - 0.01d0))
c
      sum1 = 0.d0
      sum2 = 0.d0
      sum3 = 0.d0
      sum4 = 0.d0
      do j = 1, 7
         sum1 = sum1 + ati(j)/dble(j)*(xz1**j)
         sum2 = sum2 + ati(j)/dble(j)*(xz2**j)
         sum3 = sum3 + ati(j)/dble(j+3)*(xz1**(j + 3))
         sum4 = sum4 + ati(j)/dble(j+3)*(xz2**(j + 3))
      enddo
c
      styvi = (sum2 - sum1)*(1.d0 - 0.01d0*eet)
      styvi = styvi + eet*(sum4 - sum3)
      styvi = styvi*htree*dbh/f13
c
      return
      end
******************************************************************** stem


