c   part of frt distribution
      subroutine hetk8sA
     & (lelli, ncl, shl, nzs, nthetas, thetas,
     & stdns, htr, hc1, hc2, rcr, dbh, glmp, ulg, uuu, gaps)
c     Calculation of gap probabilities in crown layer 
c        basically, a wrapper for spooi
c     PART A, unidirectional gaps for upward directions
c     PARTB, bi-directional gaps further down in the file
c     Tõenäosuste arvutamine               A. Kuusk   1.12.2000
c     Modified by Matti Mõttus (2020)
c   ------------------------------------------------- input parameters
c   lelli (logical): whether crown shape is ellipsoid
c   ncl: no. of tree classes
c   shl: shoot length
c   nzs: number of crown layers in numerical integration
c   nthetas: no. of zenith angles used (length of thetas)
c   thetas: array of zenith angles
c   stdns: stand density [1/m]
c   htr: tree height [m]
c   hc1: crown length, ell / con [m]
c   hc2: crown length, cylinder [m]
c   rcr: crown radius [m]
c   dbh: trunk diameter at brest height [m]
c   glmp: Fischer's grouping index
c   ulg: uuu / 2 (=uuu*G?)???
c   uuu: effective plant area density (leaves corrected for clumping + branches)
c   ------------------------------------------------ output parameters
c   gaps: gap probability
      implicit none
      include 'frtpar.h'
c     input
      logical lelli(nclmax)
      integer ncl, nzs, nthetas
      double precision shl(nclmax), thetas(nthetas), stdns(nclmax)
     & , htr(nclmax), hc1(nclmax), hc2(nclmax), rcr(nclmax)
     & , dbh(nclmax), glmp(nclmax), ulg(nclmax), uuu(nclmax)
c     output
      double precision gaps(nthetas)
c     local variables
      double precision rlls1, rllv1, rllv3, theta, stheta,
     & ctheta, sphi2, cphi2, calph, chs1, chs3, aasi, poodi
      double precision x1, y1, z1, x2, y2, z2, x3, y3, z3
      integer itheta
      save x1, y1, z1, x2, y2, z2
c     WWW unclear, why save is necessary, apparently these variables are never assigned to
c         it would be best to hard-code these values
      data x1, y1, z1, x2, y2, z2, x3, y3, z3 /9*0.d0/

c                                   gap probability on ground
c                                              gaps(maapinnal)
c     WWW unclear, why these values are not given directly to spooi.
      rlls1 = 0.d0
      rllv1 = 0.d0
      rllv3 = 0.d0
      chs3   = 0.d0
      sphi2  = 0.d0
      cphi2  = 1.d0
      chs1   = 0.d0
      calph  = 1.d0
      do itheta = 1, nthetas
         theta  = thetas(itheta)
         stheta = sin(theta)
         ctheta = cos(theta)
c
c    calculate aas: gap probability, Sun direction
c    WWW actually, this loops over ALL angles, including the view angle thknot(1)
c    and the following cubature (2,nknot). Transmittance over cubature is used to 
c    calculate optical effective LAI in funcf.
         call spooi
     & (lelli, ncl, ulg, shl, stdns, htr, hc1, hc2, rcr, dbh,
     & glmp, x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
     & stheta, ctheta, stheta, ctheta, sphi2, cphi2, calph, chs1, chs3,
     & aasi, poodi)
         gaps(itheta) = aasi
      enddo
      return
      end    


c --------------------------------------------------------------------
      subroutine hetk8sB
     & (lelli, ncl, shl, nzs, nthetav, thetav, thetasun, phi,
     & stdns, htr, hc1, hc2, rcr, dbh, glmp, ulg, uuu,
     & psgvu, bdgfuT, bdgfdT, btr1ukT, btr1dkT)
c     Calculation of gap probabilities in crown layer
c     Part B, bidirectional gaps
c     Tõenäosuste arvutamine               A. Kuusk   1.12.2000
c     Modified by Matti Mõttus (2020)

c   ------------------------------------------------- input parameters
c   lelli (logical): whether crown shape is ellipsoid
c   ncl: no. of tree classes
c   shl: shoot length
c   nzs: number of crown layers in numerical integration
c   nthetav: no. of zenith angles used 
c   thetav: vector of view zenith angles [rad], e.g. knots of G-L quadrature (zenith angle)
c   thetasun: sun zenith angle [rad]
c   phi: relative view azimuth angle [rad]. All directions are assumed to have the same phi
c   stdns: stand density [1/m]
c   htr: tree height [m]
c   hc1: crown length, ell / con [m]
c   hc2: crown length, cylinder [m]
c   rcr: crown radius [m]
c   dbh: trunk diameter at brest height [m]
c   glmp: Fischer's grouping index
c   ulg: uuu / 2 (=uuu*G?)???
c   uuu: effective plant area density (leaves corrected for clumping + branches)
c   ------------------------------------------------ output parameters
c   Note: compared with original frt arrays, the output arrays were transposed in 2020
c      and the letter T was added to their names
c   psgvu: bidirectional (sun, view) gap probability, upper hemisphere
c   ===WARNING== the contents from the next variables were guessed from
c     comment lines without going through the code line by line!
c   bdgfuT: Integral of bidirectional gap probability over a crown, upward direction
c   bdgfdT: Integral of bidirectional gap probability over a crown, downward direction
c   btr1ukT: probability of seeing sunlit trunk from above
c   btr1dkT: probability of seeing sunlit trunk from below

      implicit none
      include 'frtpar.h'
c  input
      logical lelli(nclmax)
      integer ncl, nzs, nthetav
      double precision shl(nclmax), thetav(nknotm), thetasun, phi,
     & stdns(nclmax), htr(nclmax), hc1(nclmax), hc2(nclmax),
     & rcr(nclmax), dbh(nclmax), glmp(nclmax), ulg(nclmax),
     & uuu(nclmax)
c  output
      double precision psgvu(nthetav), bdgfuT(nclmax,nthetav),
     & btr1ukT(nclmax,nthetav), bdgfdT(nclmax,nthetav),
     & btr1dkT(nclmax,nthetav)
c  local variables
      integer itheta, icl
      double precision
     & cphi, theta1, stheta1, ctheta1, theta2, stheta2, ctheta2,
     & rlls1, rllv1, rllv3, chs1, chs3, calph, pooui, poodi, ul,
     & sphi
      double precision x1, y1, z1, x2, y2, z2, x3, y3, z3
      save x1, y1, z1, x2, y2, z2
c     WWW unclear, why save is necessary, apparently these variables are never assigned to
c         it would be best to hard-code these values
      data x1, y1, z1, x2, y2, z2, x3, y3, z3 /9*0.d0/
c
      rlls1 = 0.d0
      rllv1 = 0.d0
      rllv3 = 0.d0
      chs1   = 0.d0
      chs3   = 0.d0
c   rlls1, rllv1, rllv3, chs1 and chs3 are between-crown hotspot parameters in spooi.f
c
c    BIDIRECTIONAL GAPS
c
      sphi   = sin(phi)
      cphi   = cos(phi)
      theta1  = thetasun
      stheta1 = sin(theta1)
      ctheta1 = cos(theta1)
      do itheta = 1, nthetav
        theta2    = thetav(itheta)
c                                p(sunlit, ground) and pooui(ground)
c                              p(sunlit, ground) ja pooui(maapinnal)
        stheta2 = sin(theta2)
        ctheta2 = cos(theta2)
        calph  = stheta2*stheta1*cphi + ctheta2*ctheta1
c      calph=cos(alpha), alpha is the angle between Sun & view directions
c
        call spooi
     & (lelli, ncl, ulg, shl, stdns, htr, hc1, hc2, rcr, dbh,
     & glmp, x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
     & stheta1, ctheta1, stheta2, ctheta2, sphi, cphi, calph,chs1,chs3,
     & pooui, poodi)
c
        psgvu(itheta) = pooui
c
        do icl = 1, ncl
c         subroutine enel fills only one element in bdgfuT, bdgfdT, btr1ukT and btr1dkT at a time:
c           that element is indicated by icl, itheta.
c         WWW it would make sense to move the loops into enel, 
c           or at least make it accept one element at a time
          ul = uuu(icl)
          call enel
     & (lelli, icl, itheta, ncl, ul, shl, stdns, htr, dbh, hc1,hc2,rcr,
     & ulg, glmp, thetasun, theta2, phi, nzs,
     & bdgfuT, bdgfdT, btr1ukT, btr1dkT)
        enddo
      enddo ! do itheta = 1, nthetav
c
      return
      end
*
************************************************************************
*
