!   part of the FRT distribution
c   this file contains spooi & friends: pi11d, pi11u, pi22d, pi22u
c     depends on rmsub.f for scone, stem
c   compile: f2py.exe -c --compiler=mingw32 -m spooi rmsub.o spooi.f

      subroutine spooi
     & (lelli, ncl, ulg, shl, stdns, htr, hc1, hc2, rcr, dbh,
     & glmp, x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
     & sthets, cthets, sthetv, cthetv, sphi, cphi, calph, chs1, chs3,
     & pooui, poodi)
!f2py intent(in) lelli, ncl, ulg, shl, stdns, htr, hc1, hc2, rcr, dbh, glmp, x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3, sthets, cthets, sthetv, cthetv, sphi, cphi, calph, chs1, chs3
!f2py intent(out) pooui, poodi
c   bidirectional probability of between-crown gaps:
c   ------------------------------------------ some input variables:
c   Sun and view directions given by the pairs [sin(t),cos(t)]:
c    [sthets, cthets] & [sthetv, cthetv]
c   if  [sthets, cthets] = [sthetv, cthetv], calculate monodirectional
c    transmittance (one direction only, no hotspot, etc.)
c   calph=cos(alpha), alpha is the angle between Sun & view directions
c
c   ----------------------------------------------output variables:
c   pooui: bidirectional gap probability, upper hemisphere
c   poodi: bidirectional gap probability, lower hemisphere
c   (calculated using Eq. (5) in manual ver. 05.2005)
c
c   scomm oli väikese sl8 korral vale.                    25.04.2002
c   pooui - ülemine, poodi - alumine poolsfäär,  A. Kuusk 22.11.2000
c   Kahe suuna avatuse to~ena"osus väljaspool võra.   A. Kuusk 22.05.1998
c   Kroonide projektsioonid asendatakse sama pindalaga ringidega.
c   (xi, yi, zi) on kiirte päikese ning vaatleja suunas üles ja alla
c   võrast väljumise punktid: i = 1 - vaatleja suunas üles,
c   i = 2  - päike, i = 3 - vaatleja suunas alla
c    ground: z=0
c
      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c     implicit double precision (a-h, o-z)
c
      logical lellib, lthsv, lelli
      dimension lelli(1:nclmax), shl(1:nclmax), stdns(1:nclmax)
      dimension ulg(1:nclmax), rcr(1:nclmax), htr(1:nclmax),
     & dbh(1:nclmax), hc1(1:nclmax), hc2(1:nclmax), glmp(1:nclmax)
c
      data eps1/.1d-3/, eps2/.2d-5/, eps3/.1d-6/
*     data pi/3.14159265358979d0/
c
      rintegr(xxx) = (1.d0 - exp(-xxx))/xxx
c
      lthsv = (abs(sthetv - sthets) .lt. eps2) .and.
     & (abs(cphi - 1.d0) .lt. eps2)
c   lthsv: whether two directions coincide
      sumu  = 0.d0
      sumd  = 0.d0
               azs   = 0.d0
      tgths = sthets/cthets
      tgthv = sthetv/cthetv
      do 200 in = 1, ncl ! end of loop at the end of the file
         rlls2  = 0.d0
         rllv2u = 0.d0
         rllv2d = 0.d0
         x8     = 0.d0
         y8     = 0.d0
         lellib = lelli(in)
         ulgi   = ulg(in)
         rcri   = rcr(in)
         htree  = htr(in)
         dbi    = dbh(in)
         hc1i   = hc1(in)
         hc2i   = hc2(in)
         shli   = shl(in)
         if (lellib) then
            cellb = hc1i*.5d0
            hbase = htree - hc1i
         else
            hbase = htree - (hc1i + hc2i)
         endif
c                             crown projection in the view direction (up)
c                                    võra projektsioon vaatesuunas (üles)
         if (z1 .lt. htree) then
            if (lellib) then
               call  pi11u (tgthv, z1, htree, cellb, rcri, szu1, vzui)
c   szu1: projected area
c   vzui: volume of the upper part of the ellipsoid
            else
               call  pi22u
     &               (tgthv, z1, htree, hc1i, hc2i, rcri, szu1, vzui)
            endif
c
            if (szu1 .gt. 0.d0) then
               rllv2u = vzui/szu1/cthetv
               azvu  = exp(-ulgi*rllv2u) ! crown transmittance (down)
c                     or the gap probability inside the crown
               sl1   = ((htree + max(z1, hbase))*.5d0 - z1)
     &          *tgthv ! height of crown center relative to z1
            else
               azvu  = 0.d0
               sl1   = 0.d0
            endif
         else
            azvu  = 0.d0
            sl1   = 0.d0
            szu1  = 0.d0
         endif
c                           crown projection in the view direction (down)
c                                    võra projektsioon vaatesuunas (alla)
         if (z3 .gt. hbase) then
            if (lellib) then
               call  pi11d (tgthv, z3, htree, cellb, rcri, szd3, vzdi)
            else
               call  pi22d
     &               (tgthv, z3, htree, hc1i, hc2i, rcri, szd3, vzdi)
            endif
c
            if (szd3 .gt. 0.d0) then
               rllv2d = vzdi/szd3/cthetv
               azvd  = exp(-ulgi*rllv2d)
               sl3  = ((hbase + z3)*.5d0 - hbase)
     &                *tgthv
            else
               azvd  = 0.d0
               sl3   = 0.d0
            endif
         else
            azvd  = 0.d0
            sl3   = 0.d0
            szd3  = 0.d0
         endif
c                             trunk projection in the view direction (up)
c                                   tu"ve projektsioon vaatesuunas (üles)
         slty1 = 0.d0
         slty2 = 0.d0
         if (z1 .lt. hbase) then
            slty1 = (hbase - z1)*tgthv
            call stem(z1, hbase, dbi, htree, styx)
            sty1  = styx*tgthv
*           sty1  = dbi*slty1
         else
            sty1  = 0.d0
         endif
c                           trunk projection in the view direction (down)
c                                   tu"ve projektsioon vaatesuunas (alla)
         slty3 = min(z3, hbase)
         zx1 = 0.d0
         zx2 = slty3
         call stem(zx1, zx2, dbi, htree, styx)
         sty3  = styx*tgthv
*        sty3  = slty3*dbi
c                                      pooi üles binoomvalemiga, 28.05.1998
c                                      esimene kordaja
         vaheg   = 1.d0 - glmp(in)
         if (abs(vaheg) .gt. eps1) then
            b11i = -log(1.d0 - (1.d0 -azvu)*vaheg)/vaheg
            ! coefficient b1j in Eq. (7) of manual (ver. 05.2005)
         else
            b11i = 1.d0 - azvu
         endif
c
         if (lthsv) then
c  the two directions are the same, no need to calculate anything for the
c   Sun direction
            b2iz  = b11i
            b12i  = b11i
            szu2  = szu1
            scomm = szu1
            sty2  = sty1
            scty  = sty1
         else
c                               crown projection in the Sun direction
c                                    võra projektsioon Päikese suunas
            if (z2 .lt. htree) then
               if (lellib) then
                  call  pi11u
     &             (tgths, z2, htree, cellb, rcri, szu2, vzui)
               else
                  call  pi22u
     &              (tgths, z2, htree, hc1i, hc2i, rcri, szu2, vzui)
               endif
               if (szu2 .gt. 0.d0) then
                  rlls2 = vzui/szu2/cthets
                  azs   = exp(-ulgi*rlls2)
                  sl2   = ((htree + max(z2, hbase))*.5d0 - z2)
     &               *tgths
               else
                  azs   = 0.d0
                  sl2   = 0.d0
               endif
            else
               azs  = 0.d0
               sl2  = 0.d0
               szu2 = 0.d0
            endif
c                               trunk projection in the Sun direction
c                                   tu"ve projektsioon thets suunas
            if (z2 .lt. hbase) then
               slty2 = (hbase - z2)*tgths
               call stem(z2, hbase, dbi, htree, styx)
               sty2  = styx*tgths
*              sty2  = slty2*dbi
            else
               sty2  = 0.d0
            endif
c                                      distance between projection centers
c                                   projektsioonitsentrite vaheline kaugus
            x8  = x2 + sl2
            y8  = y2
            x7  = x1 + sl1*cphi
            y7  = y1 + sl1*sphi
            sl8 = sqrt((x7 - x8)**2 + (y7 - y8)**2)
c                                       common part of trunk projections
c                                         tu"ve projektsioonide u"hisosa
c
            if (sty1 .gt. 0.d0 .and. sty2 .gt. 0.d0) then
               if (cphi .lt. -eps3) then
                  scty = 0.d0
               else if (abs(y1 - y2) .lt. eps1) then
                  x5  = x2 + slty2
                  x6  = x1 + slty1*cphi
                  if (abs(x1 - x2) .lt. eps1) then
                     tstphi = dbi/min(slty1, slty2)
                     if (sphi .gt. tstphi) then
                        scty = dbi**2*.5d0/sphi
                     else if ((min(x5, x6) - max(x1, x2)) .gt. 0.d0)
     &                  then
                        scty = (min(x5, x6) - max(x1, x2))*dbi
                     else
                        scty = 0.d0
                     endif
                  else
                     if (sphi .gt. eps3) then
                        scty = 0.d0
                     else if ((min(x5, x6) - max(x1, x2)) .gt. 0.d0)
     &                  then
                        scty = (min(x5, x6) - max(x1, x2))*dbi
                     else
                        scty = 0.d0
                     endif
                  endif
               else
                  scty = 0.d0
               endif
            else
               scty = 0.d0
            endif
c                                       common part of crown projections
c                                 vo~ra projektsioonide u"hisosa pindala
            ss1 = max(szu1, szu2)
            ss2 = min(szu1, szu2)
            rc1 = sqrt(ss1/pi)
            rc2 = sqrt(ss2/pi)
            if (sl8 .ge. (rc1 + rc2)) then
               scomm = 0.d0
            else if ((sl8 .le. (rc1 - rc2)) .or. sl8 .le. eps2) then
               scomm = ss2
            else
               ppp   = (rc1 + rc2 + sl8)*.5d0
               angl4 = acos((rc1**2 + sl8**2 - rc2**2)/(2.d0*rc1*sl8))
               angl8 = acos((rc1**2 + rc2**2 - sl8**2)/(2.d0*rc1*rc2))
               angl5 = (pi - angl4 - angl8)
               ss4   = angl4*rc1**2
               ss5   = angl5*rc2**2
               ss6   = sqrt(ppp*(ppp - rc1)*(ppp - rc2)*(ppp - sl8))
               scomm = ss4 + ss5 - 2.d0*ss6
            endif
c                                       between-crown hot-spot correction
c                                        ! va"line hot-spoti korrektsioon
            rlls = rlls1 + rlls2
            rllv = rllv1 + rllv2u
            rl12 = (rllv - rlls)**2 + (1.d0 - calph)*2.d0*rlls*rllv
            if (rl12 .lt. 0.d0) then
               crr2 = 0.d0
            else
               crr2 = sqrt(rl12)/shli
            endif
            if (crr2 .lt. eps1) then
               xtmp = 1.d0 - crr2*.5d0
            else
               xtmp = rintegr(crr2)
            endif
            chs2 = exp((sqrt(rllv*rlls)*xtmp - chs1)*ulgi)
c
            azvs  = azvu*azs*chs2
            if (abs(vaheg) .gt. eps1) then
               b12i = -log(1.d0 - (1.d0 - azs)*vaheg)/vaheg
               b2iz = -log(1.d0 - (1.d0 - azvs)*vaheg)/vaheg
            else
               b12i   = 1.d0 - azs
               b2iz   = 1.d0 - azvs
            endif
         endif ! if (lthsv) then
         styx = min(sty1, sty2)
         styx = min(styx, scty)
         styx = sty1 + sty2 - styx
         sumu = sumu + stdns(in)*(b11i*(szu1 - scomm) +
     &       b12i*(szu2 - scomm) + b2iz*scomm + styx)
c
c *****************
c                                      pooi alla binoomvalemiga, 28.05.1998
c                                      esimene kordaja
         if (abs(vaheg) .gt. eps1) then
            b11i = -log(1.d0 - (1.d0 -azvd)*vaheg)/vaheg
         else
            b11i = 1.d0 - azvd
         endif
c
         x7  = x3 - sl3*cphi
         y7  = y3 - sl3*sphi
         sl8 = sqrt((x7 - x8)**2 + (y7 - y8)**2)
c
c                                         tu"ve projektsioonide u"hisosa
         if (sty3 .gt. 0.d0 .and. sty2 .gt. 0.d0) then
            if (-cphi .lt. -eps3) then
               scty = 0.d0
            else if (abs(y3 - y2) .lt. eps1) then
               x5  = x2 + slty2
               x6  = x3 - slty3*cphi
               if (abs(x3 - x2) .lt. eps1) then
                  tstphi = dbi/min(slty3, slty2)
                  if (sphi .gt. tstphi) then
                     scty = dbi**2*.5d0/sphi
                  else if ((min(x5, x6) - max(x3, x2)) .gt. 0.d0) then
                     scty = (min(x5, x6) - max(x3, x2))*dbi
                  else
                     scty = 0.d0
                  endif
               else
                  if (sphi .gt. eps3) then
                     scty = 0.d0
                  else if ((min(x5, x6) - max(x3, x2)) .gt. 0.d0) then
                     scty = (min(x5, x6) - max(x3, x2))*dbi
                  else
                     scty = 0.d0
                  endif
               endif
            else
               scty = 0.d0
            endif
         else
            scty = 0.d0
         endif
c
c                                 vo~ra projektsioonide u"hisosa pindala
         ss1 = max(szd3, szu2)
         ss2 = min(szd3, szu2)
         rc1 = sqrt(ss1/pi)
         rc2 = sqrt(ss2/pi)
         if (sl8 .ge. (rc1 + rc2)) then
            scomm = 0.d0
         else if ((sl8 .le. (rc1 - rc2)) .or. sl8 .le. eps2) then
            scomm = ss2
         else
            ppp   = (rc1 + rc2 + sl8)*.5d0
            s34   = ppp*(ppp - rc1)*(ppp - rc2)*(ppp - sl8)
            sl4   = 2.d0/sl8*sqrt(s34)
            angl1 = 2.d0*asin(sl4/rc1)
            angl2 = 2.d0*asin(sl4/rc2)
            ss4   = rc1**2*.5d0*(angl1 - sin(angl1))
            ss5   = rc2**2*.5d0*(angl2 - sin(angl2))
            scomm = ss4 + ss5
         endif
c
c                                        ! va"line hot-spoti korrektsioon
         rlls = rlls1 + rlls2
         rllv = rllv3 + rllv2d
         rl12 = (rllv - rlls)**2 + (1.d0 + calph)*2.d0*rlls*rllv
         if (rl12 .lt. 0.d0) then
            crr2 = 0.d0
         else
            crr2 = sqrt(rl12)/shli
         endif
         if (crr2 .lt. eps1) then
            xtmp = 1.d0 - crr2*.5d0
         else
            xtmp = rintegr(crr2)
         endif
         chs2 = exp((sqrt(rllv*rlls)*xtmp - chs3)*ulgi)
c
         azvs  = azvd*azs*chs2
         if (abs(vaheg) .gt. eps1) then
            b12i = -log(1.d0 - (1.d0 - azs)*vaheg)/vaheg
            b2iz = -log(1.d0 - (1.d0 - azvs)*vaheg)/vaheg
         else
            b12i   = 1.d0 - azs
            b2iz   = 1.d0 - azvs
         endif
c
         styx = min(sty3, sty2)
         styx = min(styx, scty)
         styx = sty3 + sty2 - styx
         sumd = sumd + stdns(in)*(b11i*(szd3 - scomm) +
     &         b12i*(szu2 - scomm) + b2iz*scomm + styx)
c
200   continue
c
      pooui = exp(-sumu)
c
      poodi = exp(-sumd)
c
      return
      end
*
******************************************************************** spooi
*

      subroutine pi11d
     & (tgthx, zz12, htree, cellb, rcri, szdx, vzdi)
cf2py intent(in)tgthx, zz12, htree, cellb, rcri
cf2py intent(out) szdx, vzdi
c                                                 A. Kuusk  16.11.2000
c   Ellipsoidi alumise osa projektsioon szdx ja ruumala vzdi
c   ko~rgusel zz12 suunas thx 
c   tgthx = tan(thx), htree - puu ko~rgus,
c   cellb - ellipsi pooltelg, rcri - krooni raadius
c
      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c      
      data eps/.1d-4/
c
      hcrown = 2.d0*cellb
      hbase  = htree - hcrown
c                                                              tüvi
      if (zz12 .ge. 0.d0 .and. zz12 .le. hbase) then
         szdx = 0.d0
         vzdi = 0.d0
         return
      endif
c                                                   võrast kõrgemal
      if (zz12 .ge. htree) then
         vzdi = 4.d0*pi*rcri**2*cellb/3.d0
         szdx = sqrt(rcri**2 + (cellb*tgthx)**2)*pi*rcri
         return
      endif
c
      hcent = htree - cellb
      zx    = zz12 - hcent
      vzdi  = pi*rcri**2*(hcrown/3.d0 + zx - zx**3/3.d0/cellb**2)
!     MATTI: vzdi can become very slightly negative throwing an NaN downstream
      if ( vzdi .lt. eps ) then
         vzdi = 0.d0
      endif
c                                                   puutepunkti kõrgus
      if (tgthx .lt. eps) then
         z0 = 0.d0
      else
         z0 = cellb/sqrt(rcri/(cellb*tgthx)**2 + 1.d0) 
      endif
c
      rz = sqrt(cellb**2 - zx**2)*rcri/cellb
      if (abs(zx) .ge. z0) then
         if (zz12 .gt. hcent) then
c                                                   võra ülemine ots
            szdx = sqrt(rcri**2 + (cellb*tgthx)**2)*pi*rcri
         else
c                                                   võra ülemine ots
            szdx = pi*rz**2
         endif
      else
c                                                   puutepunktide vahel
         xyz  = sqrt(rcri**2 + (cellb*tgthx)**2)
         beta = acos(zx*tgthx/xyz)
         sel1 = rcri*xyz*beta - rz*zx*tgthx
         sel3 = rcri*xyz*(pi - beta) + rz*zx*tgthx
         szdx = sel1 + sel3
      endif
      return
      end
*
******************************************************************** pi11d
      subroutine pi11u
     & (tgthx, zz12, htree, cellb, rcri, szux, vzui)
cf2py intent(in) tgthx, zz12, htree, cellb, rcri
cf2py intent(out) szux, vzui
c   Projection of an ellipsoid on a horizontal surface
c   -------------- input parameters
c   tgthx: tan(t), where t = zenith angle of the ray used for projection
c   zz12: height of the surface
c   htree: tree height
c   cellb: vertical semiaxis of the ellipsoid (crown length/2)
c   rcri: horizontal semiaxis of the ellipsoid (crown radius)
c   -------------- output parameters
c   szux: projected area
c   vzui: volume of the upper part of the ellipsoid
c
c                                              A. Kuusk  16.11.2000
c   Ellipsoidi projektsioon horisontaalsele pinnale szux
c   ko~rgusel zz12 suunas thx ja u"lemise osa ruumala vzui,
c   tgthx = tan(thx), htree - puu ko~rgus,
c   cellb - ellipsi pooltelg, rcri - krooni raadius
c
      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c
      data eps/.1d-4/
c                                              zz12 is above the crown
      if (zz12 .ge. htree) then
         szux = 0.d0
         vzui = 0.d0
         return
      endif
c
      hcrown = 2.d0*cellb
      hbase  = htree - hcrown
c                                              zz12 is below the crown
      if (zz12 .ge. 0.d0 .and. zz12 .le. hbase) then
         ! return total projected area and total crown volume
         szux = sqrt(rcri**2 + (cellb*tgthx)**2)*pi*rcri
         vzui = 4.d0*pi*rcri**2*cellb/3.d0
         return
      endif
c
      hcent = htree - cellb ! height of the centre of the ellipsoid
      zx    = zz12 - hcent ! height relative to crown center
      vzui  = pi*rcri**2*(hcrown/3.d0 - zx + zx**3/3.d0/cellb**2)
c MATTI WWW: vzui could occasionally become slightly negative causing problems downstream
      if (vzui .lt. eps) then
          vzui=0.d0
      endif
        ! volume of crown above zz12
c                                                   puutepunkti kõrgus
c     height above crown center at which the ray touches the ellipsoid
      if (tgthx .lt. eps) then
         z0 = 0.d0
      else
         z0 = cellb/sqrt(rcri/(cellb*tgthx)**2 + 1.d0)
      endif
c
      rz = sqrt(cellb**2 - zx**2)*rcri/cellb
      if (abs(zx) .ge. z0) then
c                                                   võra ülemine ots
         if (zz12 .gt. hcent) then
            szux = pi*rz**2
         else
c                                                   võra alumine ots
            szux = sqrt(rcri**2 + (cellb*tgthx)**2)*pi*rcri
         endif
      else
c                                                   puutepunktide vahel
         xyz  = sqrt(rcri**2 + (cellb*tgthx)**2)
         beta = acos(zx*tgthx/xyz)
         sel1 = rcri*xyz*beta - rz*zx*tgthx
         sel2 = pi*rz**2/2.d0
         szux = sel1 + sel2
      endif
      return
      end
*
******************************************************************** pi11u
*

      subroutine pi22d
     & (tgthx, zz12, htree, hc1i, hc2i, rcri, szdx, vzdi)
cf2py intent(in) tgthx, zz12, htree, hc1i, hc2i, rcri
cf2py intent(out) szdx, vzdi
c                                       A. Kuusk   14.11.2000
c   Tasandist zz12 madalamale jääva võra osa
c   ruumala ja projektsioon suunas thx, koonus + silinder
c   tgthx = tan(thx), htree - puu ko~rgus,
c   hc1i, hc2i - koonuse ja silindri ko~rgus, rcri - krooni raadius
c   szdi - alumise osa projektsiooni pindala, vzdi - alumise osa ruumala
c
      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c      
*     data pi/3.14159265358979d0/
c
*     if (zz12 .lt. 0.d0) then
*        print *, 'Error: z = ', zz12
*        stop
*     endif
c
      hbase  = htree - (hc1i + hc2i)
c                                                   tüvi
      if (zz12 .le. hbase) then
         szdx = 0.d0
         vzdi = 0.d0
         return
      endif
c
      vcrown = rcri**2*pi*(hc1i/3.d0 + hc2i)
      scyl   = 2.d0*hc2i*rcri*tgthx
c                                                   z > htree
      if (zz12 .ge. htree) then
         vzdi   = vcrown
         call scone(hc1i, hc1i, rcri, tgthx, scon)
         scrown = scyl + scon
         szdx   = scrown
         return
      endif
c                                                   silinder
      hcbase = htree - hc1i
      if (zz12 .lt. hcbase) then
         vzdi  = rcri**2*pi*(zz12 - hbase)
         scyld = 2.d0*rcri*tgthx*(zz12 - hbase)
         szdx  = scyld + rcri**2*pi
         return
      endif
c                                                   koonus
      hcz  = htree - zz12
      rcz  = hcz*rcri/hc1i
      vzui = rcz**2*pi*hcz/3.d0
      vzdi = vcrown - vzui
      call scone(hc1i, hcz, rcri, tgthx, scdx)
      szdx = scdx + scyl
c
      return
      end
*
******************************************************************** pi22d
*
      subroutine pi22u
     & (tgthx, zz12, htree, hc1i, hc2i, rcri, szux, vzui)
!f2py intent(in) tgthx, zz12, htree, hc1i, hc2i, rcri
!f2py intent(out) szux, vzui
c                                       A. Kuusk   14.11.2000
c   Tasandist zz12 kõrgemale jääva võra osa
c   ruumala ja projektsioon suunas thx, koonus + silinder
c   tgthx = tan(thx), htree - puu ko~rgus,
c   hc1i, hc2i - koonuse ja silindri ko~rgus, rcri - krooni raadius
c   szui - ülemise osa projektsiooni pindala, vzui - ülemise osa ruumala
c
      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c      
*     data pi/3.14159265358979d0/
c
*     if (zz12 .lt. 0.d0) then
*        print *, 'Error: z = ', zz12
*        stop
*     endif
c                                                   z > htree
      if (zz12 .ge. htree) then
         szux = 0.d0
         vzui = 0.d0
         return
      endif
c
      vcrown = rcri**2*pi*(hc1i/3.d0 + hc2i)
      call scone(hc1i, hc1i, rcri, tgthx, scon)
c                                                   tüvi
      if (zz12 .lt. (htree - (hc1i + hc2i))) then
         vzui = vcrown
         szux = 2.d0*hc2i*rcri*tgthx + scon
         return
      endif
c                                                   silinder
      hcbase = htree - hc1i
      if (zz12 .lt. hcbase) then
         vzui  = rcri**2*pi*(hc1i/3.d0 + hcbase - zz12)
         scylu = 2.d0*rcri*tgthx*(hcbase - zz12)
         szux  = scylu + scon
         return
      endif
c                                                   koonus
      hcz  = htree - zz12
      rcz  = hcz*rcri/hc1i
      vzui = rcz**2*pi*hcz/3.d0
      call scone(hcz, hcz, rcz, tgthx, szux)
c
      return
      end
*
******************************************************************** pi22u
*


