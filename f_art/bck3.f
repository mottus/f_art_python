!    part of the frt distribution
!      depends on rmsub.f (rcone, rlips) and spooi.f
      subroutine bck3
     & (lelli, icl, ncl, ul, shl, stdns, htr, dbh, hc1, hc2, rcr,
     & ulg, glmp, sthets, cthets, sthetv, cthetv, sphi,
     & cphi, calph, xi, yj, zk, sxui, sxdi)
!f2py intent(in) lelli, icl, ncl, ul, shl, stdns, htr, dbh, hc1, hc2, rcr, ulg, glmp, sthets, cthets, sthetv, cthetv, sphi, cphi, calph, xi, yj, zk
!f2py intent(out) sxui, sxdi
c  compile: f2py.exe -c --compiler=mingw32 -m bck3 rmsub.o spooi.o bck3.f
c
c     1. ja"rku hajumise heleduskoefitsient (vo~rad)
c     Jo~eveer, (mai 88) Andrese progr.-st ilips xi, yj, zk
c
c    Bidirectional gap probability in a tree crown for the spherical
c    orientation of leaves. A. Kuusk, The hot spot effect in plant canopy
c    reflectance. - R. Myneny and J. Ross (Eds.), Photon-Vegetation
c    Interactions. Applications in Optical Remote sensing and Plant
c    Ecology. Springer, Berlin, 1991, 139-159.   A. Kuusk  4.01.1986
c    Muudetud 11.05.1998
c    lellips = .true. - ellipsoid, = .false. -- koonus + silinder
c    calph = cos(alph), alph - nurk kahe kiire vahel
c    shl - hot spot parameter, lehe (kasvu) suurus (m)
c    sthet? = sin(thet?), cthet? = cos(thet?)
c    sphi = sin(phi), cphi = cos(phi), phi on nurk kahte kiirt
c    la"bivate vertikaaltasandite vahel.
c    aell - võra raadius, vhelko - võra ko~rgus, htree - puu ko~rgus
c
      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c      
      logical lellib, lelli
      dimension lelli(1:nclmax)
c
      dimension shl(1:nclmax), stdns(1:nclmax), htr(1:nclmax), 
     & dbh(1:nclmax), hc1(1:nclmax), hc2(1:nclmax), rcr(1:nclmax), 
     & ulg(1:nclmax), glmp(1:nclmax)
c
      rintegr(xxx) = (1.d0 - exp(-xxx))/xxx
c
      vhelko = hc1(icl)
      lellib = lelli(icl)
      aell   = rcr(icl)
      htree  = htr(icl)
      shli   = shl(icl)
c
c  sthedv = sthetv, cthedv = -cthetv, sphid = +sphi, 
c  cphid = -cphi, calphd = -calph
c
      if (lellib) then
         cell   = .5d0*vhelko
c
         call rlips
     &    (xi, yj, zk, sthetv, cthetv, sphi, cphi, aell, cell, rlout)
         rllvu  = rlout
         call rlips
     &  (xi, yj, zk, sthetv, -cthetv, sphi, -cphi, aell, cell, rlout)
         rllvd  = rlout
c
         call rlips
     &    (xi, yj, zk, sthets, cthets, 0.d0, 1.d0, aell, cell, rlout)
         rlls = rlout
c
      else
         cell   = vhelko
         vhcyl  = hc2(icl)
c
         call rcone
     &   (xi, yj, zk, sthetv, cthetv, sphi, cphi, 
     &    aell, vhelko, vhcyl, rlout)
         rllvu  = rlout
         call rcone
     &   (xi, yj, zk, sthetv, -cthetv, sphi, -cphi, 
     &    aell, vhelko, vhcyl, rlout)
         rllvd  = rlout
c
         call rcone
     &   (xi, yj, zk, sthets, cthets, 0.d0, 1.d0, 
     &    aell, vhelko, vhcyl, rlout)
         rlls = rlout
c
      endif
c
c                                               ! hot-spoti korrektsioon
      rl12u = (rllvu - rlls)**2 + (1.d0 - calph)*2.d0*rlls*rllvu
      if (rl12u .lt. 0.d0) then
         crr0 = 0.d0
      else
         crr0 = sqrt(rl12u)/shli
      endif
      if (crr0 .lt. .001d0) then
         xtmp = 1.d0 - crr0*.5d0
      else
         xtmp = rintegr(crr0)
      endif
      chs1 = sqrt(rllvu*rlls)*xtmp
      xx  = (chs1 - rlls - rllvu)*ul*.5d0
      bgpu = exp(xx)
c
      rl12d = (rllvd - rlls)**2 + (1.d0 + calph)*2.d0*rlls*rllvd
      if (rl12d .lt. 0.d0) then
         crr0 = 0.d0
      else
         crr0 = sqrt(rl12d)/shli
      endif
      if (crr0 .lt. .001d0) then
         xtmp = 1.d0 - crr0*.5d0
      else
         xtmp = rintegr(crr0)
      endif
      chs3 = sqrt(rllvd*rlls)*xtmp
      xx  = (chs3 - rlls - rllvd)*ul*.5d0
      bgpd = exp(xx)
c                bgp - kahe suuna ava to~ena"osus võra sees
c                      bidirectional gap probability
c
c  To~ena"osuste pooi arvutamine.  (A. Jo~eveer)
c  Vaba vaatesuuna tõenäosused väljaspool võra korraga päikese suunas
c  ja vaatesuunas nii üles kui alla, pooui ja poodi, vastavalt.
c
      x1  = xi + rllvu*sthetv*cphi
      y1  = yj + rllvu*sthetv*sphi
      z1  = zk + rllvu*cthetv + htree - cell
c                     htree - puu ko~rgus, 
c                     (htree - cell) - võra koordinaadistiku alguspunkt
      x2  = xi + rlls*sthets
      y2  = yj
      z2  = zk + rlls*cthets + htree - cell
c
      x3  = xi - rllvd*sthetv*cphi
      y3  = yj + rllvd*sthetv*sphi
      z3  = zk - rllvd*cthetv + htree - cell
c
      call spooi
     & (lelli, ncl, ulg, shl, stdns, htr, hc1, hc2, rcr, dbh, 
     & glmp, x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls, rllvu, rllvd, 
     & sthets, cthets, sthetv, cthetv, sphi, cphi, calph, chs1, chs3,
     & pooui, poodi)
c
      sxui  = bgpu*pooui
      sxdi  = bgpd*poodi

c
      return
      end
*
************************************************************************
*
