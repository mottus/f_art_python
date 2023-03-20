      subroutine hetk8o
     & (ncl, thets, thetv, phi, stdns, uuu,
     & bdgfu, bdgfd, btr1uk, btr1dk, pkhair, rnlf, rtrnk, rrs, ttt,
     & bi0u, bi0d )
     

c  Single scattering of crowns as viewed from above and below
c                           [ Eq. (3) of manual ver. 09.2002]
c  võrade ühekordse hajumise heledus alt ja ülalt vaadates
c     A. Kuusk  27.11.2000
c
c IN: ncl, thets, thetv, phi, stdns, uuu, bdgfu, bdgfd, btr1uk,
c   btr1dk, pkhair, rnlf, rtrnk, rrs, ttt
c OUT:
c   bi0u: single scattering reflectance component (tree crowns)
c   bi0d: single scattering transmittance component (tree crowns)

      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c
      dimension stdns(nclmax), uuu(nclmax)
      dimension bdgfu(nclmax), btr1uk(nclmax)
      dimension bdgfd(nclmax), btr1dk(nclmax)
      dimension rnlf(nclmax), rtrnk(nclmax)
      dimension rrs(nclmax), ttt(nclmax)
c
c     data pi/3.14159265358979d0/
      data eps/.01d0/
c
      sthetv  = sin(thetv)
      cthetv  = cos(thetv)
      cphi    = cos(phi)
      sphi    = sin(phi)
c
      cthets  = cos(thets)
      sthets  = sin(thets)
      tgths   = sthets/cthets
      calph   = sthetv*sthets*cphi + cthetv*cthets
      alpha   = acos(calph)
      aalpha  = abs(alpha)
      if (aalpha .lt. eps) alpha = 0.d0
      if (alpha .lt. 0.d0) alpha = alpha + pi
      salph   = sin(alpha)
      calp2a  = cos(alpha*.5d0)
      calp2b  = sqrt((1.d0 - calph)/2.d0)
c
      bi0u    = 0.d0
      bi0d    = 0.d0
c                                  !   Vo~rade u"hekordse hajumise heleduse
c                                  !   summeerimine u"le suurusklasside
c        ---------------------------------------------------------
c        loop over tree classes:
c           single-scattering reflectance [ Eq. (3) of manual ver. 09.2002]
      do icl = 1, ncl
            uli    = uuu(icl)
            rdlf   = rrs(icl)
            tleaf  = ttt(icl)
            rni    = rnlf(icl)
c          ------ upper hemisphere
c           calculate single scattering reflectance from leaves
c           Gamma,  spherical orientation, upper hemisphere
            gammr = (salph + (pi - alpha)*calph)/(3.d0*pi)
            gammt = (salph - alpha*calph)/(3.d0*pi)
            gammd = rdlf*gammr + tleaf*gammt
c           fresnel reflectance
            call gmfres(calp2a, rni, pkhair, rfr)
c
            gammout = gammd + rfr*.125d0
            bc1u    = bdgfu(icl)*gammout*uli*stdns(icl)/
     &                (cthets*cthetv)
c  WWW this should be Eq. (3) of manual ver. 09.2002. Since the only variable
c    depending on x,y,z is bidirectional gap probability, it is integrated
c    beforehand (bdgfu). but WHY IS IT DIVIDED BY cthetv? Is it because of
c    the definition of Gamma?

c           calculate trunk reflectance
c           Gamma, vertical orientation (tty = 0, rfr = 0, upper hemisphere)
            gammout = rtrnk(icl)*(sphi - (phi - pi)*cphi)/(2.d0*pi)
            btr1u   = btr1uk(icl)*gammout*stdns(icl)*
     &                tgths*sthetv/cthetv
c
            bi0u    = bi0u + bc1u + btr1u
c
c          ------ lower hemisphere
            gammd = rdlf*gammt + tleaf*gammr
            call gmfres(calp2b, rni, pkhair, rfr)
            gammout = gammd + rfr*.125d0
            bc1d    = bdgfd(icl)*gammout*uli*stdns(icl)/
     &                (cthets*cthetv)
c           Gamma,  spherical orientation, lower hemisphere
            gammout = rtrnk(icl)*(sphi - phi*cphi)/(2.d0*pi)
c
            btr1d   = btr1dk(icl)*gammout*stdns(icl)*
     &                tgths*sthetv/cthetv
            bi0d    = bi0d + bc1d + btr1d


      enddo
      return
      end
*
***********************************************************************
*
