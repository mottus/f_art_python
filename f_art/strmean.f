!   part of the fortran-only version of FRT
! NOTE: not used by the python version of frt, f2py definitions provided for testing
      subroutine strmean
     & (lelli, ncl, stdns, htr, hc1, hc2, rcr, dbh,
     & rmass, slwcl, rlai, rbai, clmpst, clmpsh, glmp, ulg, uuu,
     & efflai, sntr, hmtree, vhekm, vhcilm, rmcrown, dbhmean, rmassm,
     & slwm, vliit, cano, tlty, tlaief, tlai, tbai)
!f2py intent(in) lelli, ncl, stdns, htr, hc1, hc2, rcr, dbh, rmass, slwcl, rlai, rbai, clmpst, clmpsh, glmp
!f2py intent(out) ulg, uuu, efflai, sntr, hmtree, vhekm, vhcilm, rmcrown, dbhmean, rmassm, slwm, vliit, cano, tlty, tlaief, tlai, tbai
c  depends on pi11u, pi22u -- available in spooi.f, which depends on rmsub.o
c  compile with  f2py.exe -c --compiler=mingw32 -m strmean spooi.o rmsub.o strmean.f
c
c       THE MEAN VALUES OF STRUCTURE PARAMETERS   A. KUUSK  23.09.1995
c   ------------------------------------------------- input parameters
c   lelli (logical): whether crown shape is ellipsoid
c   ncl: no. of tree classes
c   stdns: stand density
c   htr: tree height
c   hc1: crown length (ellipsoid)
c   hc2: crown length, cylinder [m]
c   rcr: crown radius [m]
c   dbh: dbh [m]
c   rmass: dry leaf wght [kg/tree]
c   slwcl: specific leaf weight [gm^-2]
c   rlai: LAI
c   rbai: BAI
c   clmpst: tree distr. parameter
c   clmpsh: shoot shading coef
c   glmp: Fischer's grouping index
c   ------------------------------------------------- output parameters
c   ulg: uuu / 2 ???
c   uuu: effective plant area density (leaves corrected for clumping + branches)
c   efflai: !   LAI of a random canopy that causes
c          extinction similar to  that of the stand at thseff (= 40°)
c   sntr: total number of trees per m^2
c   hmtree: mean tree height
c   vhekm: mean crown length (ellipsoid/cone)
c   vhcilm: mean crown length (cylinder)
c   rmcrown: mean crown radius
c   dbhmean: mean dbh
c   rmassm: mean leaf mass per m^2
c   slwm: mean SLW
c   vliit: crown cover
c   cano: canopy cover
c   tlty: total unprojected below-crown trunk area per m^2 (below crown, calculated as pi*dbh*l0/2)
c   tlaief: total effective LAI (corrected for shoot-level clumping)
c   tlai: total lai
c   tbai: total BAI

      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c
      logical lelli
c
      dimension lelli(1:nclmax)
      dimension stdns(1:nclmax), htr(1:nclmax), hc1(1:nclmax),
     &  hc2(1:nclmax), rcr(1:nclmax), dbh(1:nclmax), rmass(1:nclmax),
     &  slwcl(1:nclmax), rlai(1:nclmax), rbai(1:nclmax),
     &  clmpst(1:nclmax), clmpsh(1:nclmax),
     &  ulg(1:nclmax), uuu(1:nclmax), glmp(1:nclmax), clai(1:nclmax)
      dimension v0hm(1:nclmax)
c
      data gsf/.5d0/, thseff/.698d0/
c                                    gsf - G_spherical = 0.5
*     data pi/3.14159265358979d0/
c
c                     (********** keskmiste arvutamine  **************)
c                                 structure parameters
      sntr    = 0.d0
      hmtree  = 0.d0
      vhekm   = 0.d0
      vhcilm  = 0.d0
      rmcrown = 0.d0
      dbhmean = 0.d0
      rmassm  = 0.d0
      slwm    = 0.d0
      vliit   = 0.d0
      cano    = 0.d0
      ruum    = 0.d0  ! WWW not used anywhere?
      tlaief  = 0.d0
      tlai    = 0.d0
      tbai    = 0.d0
      tlty    = 0.d0
c
      do icl = 1, ncl
         stdi    = stdns(icl)
         sntr    = sntr + stdi ! total number of trees per m^2
         ri2     = rcr(icl)**2
         vhi     = hc1(icl) + hc2(icl) ! total crown length
         hmtree  = hmtree + stdi*htr(icl) ! weighted sum of tree height
         vhekm   = vhekm + stdi*hc1(icl) ! weighted sum of crown length (elli)
         vhcilm  = vhcilm + stdi*hc2(icl) ! weighted sum of crown length (cyl)
         if (lelli(icl)) then
            voi   = 2.d0*pi/3.d0*ri2*hc1(icl) ! crown volume
         else
            voi   = pi/3.d0*ri2*(hc1(icl) + 3.d0*hc2(icl)) ! crown volume
         endif
         ruum      = ruum + voi*stdi ! total crown volume per m^2
         v0hm(icl) = voi ! crown volume
         rmcrown   = rmcrown + stdi*rcr(icl) ! weighted sum of crown radii
         dbhmean   = dbhmean + stdi*dbh(icl) ! weighted sum of stem diam.
         rmassm    = rmassm  + stdi*rmass(icl) ! leaf mass per m^2
         slwm      = slwm    + stdi*slwcl(icl) ! weightd sum of SLW
         r2stdi    = stdi*ri2
         vliit     = vliit + r2stdi ! weighted sum of crown radii^2
         cano      = cano + clmpst(icl)*r2stdi ! used for calculating canopy closu
         tlty      = tlty + stdi*dbh(icl)*.5d0*(htr(icl) - vhi)
                                   ! total projected unshaded trunk area, projection factor 1/2
         rlaief    = rlai(icl)*clmpsh(icl) ! effective LAI
         claief    = rlaief + rbai(icl) ! effective PAI
         clai(icl) = claief ! effective PAI
         tlaief    = tlaief + rlaief ! total effective LAI
         tlai      = tlai  + rlai(icl) ! total LAI
         tbai      = tbai + rbai(icl) ! total BAI
         uu        = claief/v0hm(icl)/stdns(icl) ! effective plant area density
c WWW why not uu = claief/stdi/voi ?
         ulg(icl)  = uu/2.d0
         uuu(icl)  = uu  ! effective plant area density
      enddo
c
      hmtree  = hmtree/sntr
      vhekm   = vhekm/sntr
      vhcilm  = vhcilm/sntr
      ruum    = ruum/sntr
      rmcrown = rmcrown/sntr
      dbhmean = dbhmean/sntr
      rmassm  = rmassm/sntr
      slwm    = slwm/sntr
      vliit   = pi*vliit
      cano    = 1.d0 - exp(-pi*cano)
      tlty    = pi*tlty
c
********************  Calculation of LAI(eff) (09.2002) *****************
c                     thseff = 40° = .698 rad
         thets   = thseff
         cthets  = cos(thets)
         tgths   = sin(thets)/cthets
         alsc    = 0.d0
         zz12    = 0.d0
c
         do icl = 1, ncl
            rcri   = rcr(icl)
            htree  = htr(icl)
            hc1i   = hc1(icl)
            cellb  = hc1i/2.d0
            hc2i   = hc2(icl)
            vaheg  = 1.d0 - glmp(icl)
            if (lelli(icl)) then
               call pi11u
     &          (tgths, zz12, htree, cellb, rcri, szu1, vzui)
            else
                call pi22u
     &          (tgths, zz12, htree, hc1i, hc2i, rcri, szu1, vzui)
            endif
            if (szu1 .gt. 0.d0) then
               a1gr = exp(-gsf*clai(icl)/(stdns(icl)*szu1*cthets))
c                                      transmittance of a single crown
            else
               a1gr = 1.d0
            endif
            if (abs(vaheg) .gt. 0.1d-3) then
               b10i  = -log(1.d0 - (1.d0 - a1gr)*vaheg)/vaheg
c                               parameter c in eq. (8) in Nilson (1999)
            else
               b10i  = 1.d0 - a1gr
            endif
            alsc  = alsc + stdns(icl)*szu1*b10i
         enddo
         efflai =  cthets*alsc / gsf !   LAI of a random canopy that causes
*          extinction similar to  that of the stand at thseff (= 40°)
* WWW  in the original version, the 2 (1/G) was missing
      return
      end
*
************************************************************************
*
