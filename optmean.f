      subroutine optmean
     & ( ncl, stdns, htr, hc1, hc2, dbh, rlai, rbai, clmpsh, tlty,
     & tlaief, tbai, rlfcl, tlfcl, rnlf, rbrnc, rtrnk,
     & rleff, tleff, rneff, rty, rrs, ttt, utot, rteff, tteff)
c   part of the fortran-only version of FRT
! NOTE: not used by the python version of frt, f2py definitions provided for testing
!f2py intent(in) ncl, stdns, htr, hc1, hc2, dbh, rlai, rbai, clmpsh, tlty, tlaief, tbai, rlfcl, tlfcl, rnlf, rbrnc, rtrnk 
!f2py intent(out) rleff, tleff, rneff, rty, rrs, ttt, utot, rteff, tteff
c  compile with  f2py.exe -c --compiler=mingw32 -m optmean  optmean.f
c
c  the mean values of optical parameters    A. Kuusk 23.09.1995
c
* IN:
*  ncl: no. of canopy classes
*  stdns: <- XC( , 1) stand density
*  htr: <- XC( , 2) tree height
*  hc1: <- XC( , 3) crown length, ellipse |cone
*  hc2: <- XC( , 4) crown length, cyl
*  dbh: <- XC(jcl, 6)*1.d-2
*  rlai: raw LAI of crown classes, calculated from input data in comprt.f
*  rbai: raw BAI, calculated from input data in comprt.f
*  clmpsh: <- XC( ,10) shoot clumping parameter
*  tlty: <- strmean  total unshaded trunk area per m^2
*  tlaief: <- strmean  total effective LAI (corrected for clumping)
*  tbai:  <- strmean  total BAI
*  rlfcl: diffuse leaf reflectance (calculated in comprt.f)
*  tlfcl: leaf transittance (calculated in comprt.f)
*  rnlf: corrected wax refractive index (calculated in comprt.f)
*  rbrnc:  branch reflectance
*  rtrnk:  trunk reflectance
* OUT:
*  rleff: average reflectance of leaves+branches
*  tleff: average tranmsmittance of leaves+branches
*  rneff: average wax refractive index
*  rty: average trunk reflectance
*  rrs: average diffuse reflectace of a tree class, excl. trunks
*  ttt: ???
*  utot: total effective PAI, leaves + branches + stems
*  rteff: average effective reflectance of canopy elements
*  tteff: average leaf diffuse reflectance ?

      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c
      dimension stdns(1:nclmax), htr(1:nclmax), hc1(1:nclmax),
     & hc2(1:nclmax), dbh(1:nclmax), rlai(1:nclmax), rbai(1:nclmax),
     & clmpsh(1:nclmax), rlfcl(1:nclmax), tlfcl(1:nclmax),
     & rnlf(1:nclmax), rbrnc(1:nclmax), rtrnk(1:nclmax)
      dimension ttt(nclmax), rrs(nclmax)
c
*     data pi/3.14159265358979d0/
c
c                     (********** keskmiste arvutamine  **************)
c                                 optical parameters
c  utot - kogu efektiivne LAI (lehed + oksad + tyved)
c  rteff - kogu efektiivne r_L, tteff - kogu effektiivne t_L
c
*
      sum2   = 0.d0
      do i = 1, ncl
         stdi    = stdns(i)
         htree   = htr(i)
         vhi     = hc1(i) + hc2(i)
         cli     = stdi*dbh(i)*.5d0*(htree - vhi)
         sum2    = sum2 + cli*rtrnk(i)
      enddo
*
      rty    = 0.d0
      if (tlty .ne. 0.d0) rty = sum2*pi/tlty ! average trunk reflectance
c
      sul0   = 0.d0  ! sum of rdlf
      sumn   = 0.d0
      sut    = 0.d0
      rteff  = 0.d0
      tteff  = 0.d0
c
      do i5 = 1, ncl
         rlaief       = rlai(i5)*clmpsh(i5) ! LAI corrected for clumping
         baicl        = rbai(i5) ! BAI
         claief       = rlaief + baicl
         rdlf         = (rlaief*rlfcl(i5) + baicl*rbrnc(i5))
                 ! diffuse reflectance of a tree class elements, excl. trunks
                 !   multiplied by effective area
         rrs(i5) = rdlf/claief ! average diffuse reflectace of a tree
          ! class, excl. trunks

c
         rteff        = rteff + (rlfcl(i5) +
     &     ((1.d0 - rnlf(i5))/(1.d0 + rnlf(i5)))**2)*rlaief +
     &     rbrnc(i5)*baicl
            ! sum of weighted reflectance of tree classes (dir+dif), excl. trunks
            ! includes diffuse and specular components
         ttl          = rlaief*tlfcl(i5)
         tteff        = tteff + ttl
                 ! sum of weighted leaf diffuse transmittance of tree classes
         ttt(i5) = ttl/claief ! ??? claief is effective PAI (LAIeff+BAI)
c
         sul0    = sul0 + rdlf
           ! sum of weighted diffuse reflectances (dir+dif), excl trunks
           !   multiplied by effective area
         sut     = sut + ttl
           ! sum of weighted leaf diffuse transmittance of tree classes
* WWW sut = tteff ???
         sumn    = sumn + rlaief*rnlf(i5)
           ! sum of weighted wax refractive indices
           !   multiplied by effective area
      enddo
c
      suml  = tlaief + tbai
      utot  = suml + tlty ! total effective LAI, leaves + branches + stems
      rteff = (rteff + rty*tlty)/utot ! average effective reflectance of canopy elements
*
      tteff = tteff/utot ! average diffuse leaf reflectance

      rneff = 1.d0
      if (tlaief .ne. 0.d0) rneff = sumn/tlaief ! average wax refractive index
      tleff = sut/suml ! average diffuse transmittance of canopy elements

      r12   = sul0/suml
      rleff = r12 + ((rneff - 1.d0)/(rneff + 1.d0))**2
         ! average reflectance of leaves+branches
         ! WWW: the second term on RHS should be multiplied by sum(rlaief)/suml
         !   now, specular reflectance is also added to branches
c
      return
      end
*
************************************************************************
*
