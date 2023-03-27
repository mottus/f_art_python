c  part of frt python distribution
c     depends on spooi.f, bck3f, and indirectly rmsub.f
c            (for scone, stem in spooi.f; rlips, rcone in bck3.f)
c  compile: f2py.exe -c --compiler=mingw32 -m enel3 spooi.o bck3.o rmsub.o enel3.f
      subroutine  enel
     & (lelli, icl, iths, ncl, ul, shl, stdns, htr, dbh, hc1, hc2, rcr,
     & ulg, glmp, thets, thetv, phi, nzs,
     & bdgfuT, bdgfdT, btr1ukT, btr1dkT)
!f2py intent(in) lelli, icl, iths, ncl, ul, shl, stdns, htr, dbh, hc1, hc2, rcr, ulg, glmp, thets, thetv, phi, nzs
!f2py intent(inout) bdgfuT, bdgfdT, btr1ukT, btr1dkT
c  Bidirectional gap probability inside a crown
c   ----------------output variables (just guesses :( )
c   Note: compared with original frt arrays, the output arrays were transposed in 2020
c      and the letter T was added to their names
c  bdgfuT: Integral of bidirectional gap probability over a crown, upward direction
c  bdgfdT: Integral of bidirectional gap probability over a crown, downward direction
c  btr1ukT: probability of seeing sunlit trunk from above
c  btr1dkT: probability of seeing sunlit trunk from below
c
c       vo~ralt hajuv energia                 A. Kuusk     29.08.2002
c   Integreerib kahe suuna avatuse tõenäosuse kroonil bgp*pooi,
c   3D integraal ellipsoidis arvutatakse Hammer-Stroudi kubatuuriga,
c   Hammer, P.C., Stroud, A.H. Numerical evaluation of multiple
c   integrals II. Math. Tables and other Aids Comp. 12 (1958), No 64, 272-280
c   ja koonuses ning silindris Gauss-Legendre'i kvadratuuriga z-koordinaadi
c   järgi ning Mysovskikh (1981) algoritmiga 16.33 x-y-tasandis.
c   Kvadratuuride sõlmed ja kaalud on include-failis volint.h.
c   Tõenäosused bgp ja pooi arvutatakse bck3-s, integraalid
c   salvestatakse massiivi bdgf. Valgustatud tüve nägemise
c   tõenäosused salvestatakse massiivi btr1ik.
c
      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c
      logical lelli, lellips
      dimension lelli(1:nclmax)
c
      dimension shl(1:nclmax), stdns(1:nclmax), htr(1:nclmax),
     & dbh(1:nclmax), hc1(1:nclmax), hc2(1:nclmax), rcr(1:nclmax),
     & ulg(1:nclmax), glmp(1:nclmax)
      dimension bdgfuT(nclmax,nknotm), btr1ukT(nclmax,nknotm)
      dimension bdgfdT(nclmax,nknotm), btr1dkT(nclmax,nknotm)
c
      data xi, x2, x3, yj, y2, y3, rlls1, rllv1, rllv3,
     &     chs1, chs3, cell /12*0.d0/
*     data pi2/1.570796326794895d0/
c

      sthets = sin(thets)
      cthets = cos(thets)
      sthetv = sin(thetv)
      cthetv = cos(thetv)
      sphi   = sin(phi)
      cphi   = cos(phi)
      calph  = sthetv*sthets*cphi + cthetv*cthets
c         alph = acos(calph) on nurk kiirte r1 ja r2 vahel (rv ja rs vahel)
c
      lellips = lelli(icl)
      aell    = rcr(icl)
      a2      = aell**2
      sumu    = 0.d0
      sumd    = 0.d0
      htree   = htr(icl)
      ntstz   = 2*nzs
c                                                       ellipsoid
      if (lellips) then
         hkroon = hc1(icl)
         cell   = hkroon*.5d0
         do itst = 1, netst
            xtsti = xetst(itst)*aell
            ytsti = yetst(itst)*aell
            ztsti = zetst(itst)*cell
            atsti = aetst(itst)*a2*cell
            call bck3
     & (lelli, icl, ncl, ul, shl, stdns, htr, dbh, hc1, hc2, rcr,
     & ulg, glmp, sthets, cthets, sthetv, cthetv, sphi,
     & cphi, calph, xtsti, ytsti, ztsti, sxui, sxdi)

            sumu = sumu + sxui*atsti
            sumd = sumd + sxdi*atsti
         enddo
         
         bdgfuT(icl,iths) = sumu
         bdgfdT(icl,iths) = sumd
      else
         vhelko = hc1(icl)
         vhcil  = hc2(icl)
         hkroon = vhelko + vhcil
         do itst = 1, ntstz
c                                                        cylinder
            if (vhcil .gt. 0d0) then
               sxyu   = 0.d0
               sxyd   = 0.d0
               aztsti = acztst(itst)/2.d0*vhcil
               ztsti  = (zctst(itst) - 1.d0)/2.d0*vhcil
               rcirc  = aell
               do jtst = 1, nctst
                  xtsti = xctst(jtst)*rcirc
                  ytsti = yctst(jtst)*rcirc
                  atsti = actst(jtst)*rcirc**2
                  call bck3
     & (lelli, icl, ncl, ul, shl, stdns, htr, dbh, hc1, hc2, rcr,
     & ulg, glmp, sthets, cthets, sthetv, cthetv, sphi,
     & cphi, calph, xtsti, ytsti, ztsti, sxui, sxdi)
                  sxyu = sxyu + sxui*atsti
                  sxyd = sxyd + sxdi*atsti
               enddo
               sumu = sumu + sxyu*aztsti
               sumd = sumd + sxyd*aztsti
            endif
c                                                       cone
            if (vhelko .gt. 0d0) then
               sxyu   = 0.d0
               sxyd   = 0.d0
               aztsti = acztst(itst)/2.d0*vhelko
               ztsti  = (zctst(itst) + 1.d0)/2.d0*vhelko
               rcirc  = aell*(1.d0 - ztsti/vhelko)
               do jtst = 1, nctst
                  xtsti = xctst(jtst)*rcirc
                  ytsti = yctst(jtst)*rcirc
                  atsti = actst(jtst)*rcirc**2
                  call bck3
     & (lelli, icl, ncl, ul, shl, stdns, htr, dbh, hc1, hc2, rcr,
     & ulg, glmp, sthets, cthets, sthetv, cthetv, sphi,
     & cphi, calph, xtsti, ytsti, ztsti, sxui, sxdi)
                  sxyu = sxyu + sxui*atsti
                  sxyd = sxyd + sxdi*atsti
               enddo
               sumu = sumu + sxyu*aztsti
               sumd = sumd + sxyd*aztsti
            endif
         enddo
         bdgfuT(icl,iths) = sumu
         bdgfdT(icl,iths) = sumd
      endif
c
c
c
c    Valgustatud tu"ve na"gemise to~ena"osus
c
      btr1u = 0.d0
      btr1d = 0.d0
      utyu  = 0.d0
      utyd  = 0.d0
c                       zvari - krooni varju ko~rgus
      tgths = sthets/cthets
      if (lellips) then
         zv = sqrt(a2 + (cell*tgths)**2)/tgths
      else
         zv = aell/tgths
      endif
      zvari = htree - hkroon - zv
c
      if (zvari .gt. 0.d0) then
         if (thetv .lt. thets) then
            tgthv = sthetv/cthetv
            if (lellips) then
               zv = sqrt(a2 + (cell*tgthv)**2)/tgthv
            else
               zv = aell/tgthv
            endif
            zview = max((htree - hkroon - zv), 0.d0)
         else
            zview = zvari
         endif
c
c                               ! na"htava valgustatud tüve pindala uty
         utyd = dbh(icl)*pi/2.d0*zvari
         utyu = dbh(icl)*pi/2.d0*zview
c                               ! pooi integreerimine [0, zvari]
         dzenel = hkroon/ntstz
         if (zvari .le. dzenel) then
            dztrunk = zvari
         else
            dztrunk = dzenel
         endif
         nz = int(zvari/dztrunk + .5d0)
         nview = int(zview/dztrunk + .5d0) + 1
         dztrunk = zvari/nz
         xi      = 0.d0
         x2      = xi
         x3      = xi
         yj      = 0.d0
         y2      = yj
         y3      = yj
         rlls1   = 0.d0
         rllv1   = 0.d0
         rllv3   = 0.d0
         chs1    = 0.d0
         chs3    = 0.d0
c
         if (nview .lt. nz) then
            do 40 k = 1, nz + 1
               zk   = (k - 1)*dztrunk
               z2   = zk
               z3   = zk
c
      call spooi
     & (lelli, ncl, ulg, shl, stdns, htr, hc1, hc2, rcr, dbh,
     & glmp, xi, yj, zk, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
     & sthets, cthets, sthetv, cthetv, sphi, cphi, calph, chs1, chs3,
     & pooui, poodi)
c
               if (k .eq. 1 .or. k .eq. nz) poodi = poodi/2.d0
               btr1d = btr1d + poodi
               if (k .eq. 1 .or. k .eq. nview) pooui = pooui/2.d0
               if (k .le. nview) btr1u = btr1u + pooui
40          continue
            btr1d = btr1d/nz
            btr1u = btr1u/nview
         else
            do 41 k = 1, nz + 1
               zk   = (k - 1)*dztrunk
               z2   = zk
               z3   = zk
c
      call spooi
     & (lelli, ncl, ulg, shl, stdns, htr, hc1, hc2, rcr, dbh,
     & glmp, xi, yj, zk, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
     & sthets, cthets, sthetv, cthetv, sphi, cphi, calph, chs1, chs3,
     & pooui, poodi)
c
               if (k .eq. 1 .or. k .eq. nz) then
                  poodi = poodi/2.d0
                  pooui = pooui/2.d0
               endif
               btr1d = btr1d + poodi
               btr1u = btr1u + pooui
41          continue
            btr1d = btr1d/nz
            btr1u = btr1u/nz
         endif
      endif
      btr1ukT(icl,iths) = btr1u*utyu
      btr1dkT(icl,iths) = btr1d*utyd     
c
      return
      end
*
************************************************************************
*
