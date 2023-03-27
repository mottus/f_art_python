!    part of the frt distribution
!  not used by python version of FRT, but can be used via frt_wrapper.py
	  subroutine comprt(XC, IC, SPCIN, XOUT, SPCOUT)      
!f2py intent(in) XC, IC, SPCIN
!f2py intent(inout) XOUT, SPCOUT
c  compile with sth like
c      f2py -c  --compiler=mingw32 -m comprt bck3.o bgrdd.o enel3.o hetk8o.o hetk8s.o layer.o optmean.o rmsub.o spooi.o strmean.o twostr.o comprt.f
c  --- old func.f ----
c        XC -- forest and geometry parameter array, floating-point parameters
c        IC -- forest parameter array, integer parameters
c        SPCIN -- forest element spectral properties, matrix of double precision
c result: 
c        XOUT: scalar floating-point outputs (spectrally invariant forest stand optical characteristics)
c        SPCOUT: spectrum array outputs (spectral reflectance and transmittance factors)
c
c               A. Kuusk  24.10.1995 Avignonis 
c                   + modifications later
c  Modified by Matti Mõttus 2005-2020

      implicit none
      include 'frtpar.h'
      
      double precision XC(nclmax,*), SPCIN(nspchnl,*)
      double precision XOUT(*), SPCOUT(nspchnl,*)
      integer IC(*)

c   IC elements
c  1: no. of canopy classes (old n_cl)
c  2: correction wavelength index (not used here)
c  3: no. of crown layers in numerical integration (old n_zs)
c  4: logical flag l_opt    : whether to (re)compute element mean optical properties
c  5: logical flag l_struc : whether to (re)compute canopy overall structural characteristics
c  6: logical flag l_dirv   : whether to (re)compute everything related to view direction (gaps etc.)
c  7: logical flag l_dirs   : whether to (re)compute everything related to solar direction (gaps etc.)
c  8: logical flag l_grnd   : whether to (re)load forest floor spectra
c  9: logical flag l_refl   : whether to compute reflectance factors; if false,
c         element optical and canopy structural properties will be computed, depending on other flags, and output
c  10 & 11: start and end wavelength indices, respectively. Indices start from 1 (fortran77 convention)
c  11+i to 11+nclmax: are crown shapes ellipsoids? if >0, then yes (old l_elli)
c  12+nclmax: ijob (used to set l_flux)
c  13+nclmax: nquad_t, quadrature knots over theta (polar, zenith angle)
c  14+nclmax: nquad_p, quadrature knots over phi (azimuth)

c   XOUT elements
c  1 total number of trees per m^2, sntr
c  2 mean tree height, hmtree
c  3 mean crown length (ellipsoid, or conical part), vhekm
c  4 mean crown length of cylindrical crown part, vhcilm
c  5 mean crown radius, rmcrown
c  6 mean dbh, dbhmean
c  7 mean leaf mass per m^2, rmassm
c  8 mean SLW, slwm
c  9 total BAI (branch area index) so that PAI = tlai+tbai, tbai
c 10 total effective LAI, corrected for shoot level clumping, tlaief
c 11 total lai (leaves only), tlai
c 12 crown cover (sum of crown areas / stand area), vliit
c 13 canopy cover, cano
c 14 Diffuse non-interceptance, i-i0 integrated over hemisphere DIFN
c 15 Effective (optical) PAI as inverted from Miller's equation optPAI
c 16 Recollision probability estimate by Stenberg (2007) p_Pola
c 17 gaps in view direction gaps_v(1)
c 18 gaps in solar direction gaps_s(1)
c 19 bidirectional gaps in solar direction psgvu(1)
c 20 effective lai at 40 deg view angle efflai

c  SPCOUT columns
c  1: hemispherical-directional  reflectance factor
c  2: hemispherical-directional transmittance factor
c  3: partial reflectance factor: first-order direct by crowns
c  4: partial reflectance factor: first-order direct by ground
c  5: partial reflectance factor: higher-order reflectance + first-order diffuse sky by crowns
c  6: partial reflectance factor: higher-order reflectance + first-order diffuse sky by ground
c  7: partial transmittance factor: first-order transmittance
c  8: rteff: average reflectance of leaves+branches+trunks at current wavelength
c  9: tteff: average tranmsmittance of leaves+branches+trunks at current wavelength
 
c  XC elements. Column refers to second index
c stdns: stand density [m^-2], XC column 1
c htr: tree height [m], XC column 2
c hc1: crown length, ellipse or conical part [m], XC column 3
c hc2: crown length, cylinder [m], XC column 4
c rcr: crown radius [m], XC column 5
c dbh: trunk diameter [cm], XC column 6. When read into dbh, converted to [m]
c rmass: dry leaf weight [kg/tree], XC column 7
c slwcl: specific leaf weight, SLW [gm^-2], XC column 8
c BAI/LAI ratio: XC column 9, NOT stored in a variable
c clmpst: tree distribution parameter, cB in eq. (11) in Nilson 1999, XC column 10
c clmpsh: shoot shading coefficient (ssc), XC column 11
c crncl: wax refractive index ratio, XC column 12
c shl: hot spot parameter, shoot length (m), XC column 13
c positions 14-20 unused
c thetv: view zenith angle [rad], XC(1,21)
c phiv: view azimuth angle [rad], XC(1,22)
c thets: solar zenith angle [rad], XC(1,23)
c pkhair: leaf hair optical index (vaguely defined, set to unity), XC(1,24) [was p_khair]
c age: stand age (not used in computations) XC(1,25)
c thetv increment [rad] (not used in computations) XC(1,26)
c 
c In original frt, columns from 27 described 2-layer understory soil for the whole stand and were deleted.
c    Originally, 12-24 were leaf model parameters and were deleted, elements 25-26 (crncl, shl) are now 12-13

c ------  SPCIN elements (columns). Column refers to second index
c 1  wavelength (nm) 
c 2  hemispherical-hemispherical reflectance of forest floor (=albedo) (rddgrou)
c 3  directional-hemispherical reflectance of forest floor (rsdgrou)
c 4  hemispherical-directional reflectance of forest floor (rdsgrou)
c 5  BRDF of forest floor (directional-directional reflectance) (rsogrou)
c 6  S/Q ratio (direct/total irradiance) at TOC (s_qarr)
c 7  leaf wax refractive index
c 8          to 8+nclmax-1    leaf reflectance (p_rarr) for each individual tree class
c 8+nclmax   to 8+2*nclmax-1  leaf transmittance (p_tarr)
c 8+2*nclmax to 8+3*nclmax-1  leaf abaxial reflectance [NOT USED]
c 8+3*nclmax to 8+4*nclmax-1  branch reflectance (b_rrarr) 
c 8+4*nclmax to 8+5*nclmax-1  trunk reflectance (t_rrarr)
c 8+5*nclmax   correction factor for diffuse fluxes (to correct for the error in estimating 1st order scattering) (c_fact)

c ------ input tree class-specific stand parameters
      double precision stdns(nclmax), htr(nclmax), hc1(nclmax),
     & hc2(nclmax), rcr(nclmax), dbh(nclmax), rmass(nclmax),
     & slwcl(nclmax), clmpst(nclmax), clmpsh(nclmax),
     & crncl(nclmax), shl(nclmax)     
      save stdns, htr, hc1, hc2, rcr, dbh, rmass, slwcl,
     & clmpst, clmpsh, crncl, shl
c  described above

c ------  GEOMETRY ------
      double precision thetv, phiv, thets
      double precision sqratio
      save thetv, phiv, thets
c thetv: view zenith angle, radians (was in t_hetv)      
c phiv: relative view azimuth angle [rad] (was p_hi)
c thets: solar zenith angle [rad]
c sqratio: direct to total irradiance ratio for current wavelength
     
c ------ DERIVED TREE CLASS-SPECIFIC STAND PARAMETERS
      double precision rlai(nclmax), rbai(nclmax), glmp(nclmax)
      integer ncl, nzs
      logical lelli(nclmax)
      save rlai, rbai, glmp, ncl, lelli
c rlai: LAI (no branches etc.) - from crown leaf area 
c rbai: branch area index (BAI) - from BAI/LAI ratio
c glmp: Fisher's grouping index  - tree distribution parameter
c ncl: number of canopy classes (from IC(1))
c nzs: no. of crown layers in numerical integration (old n_zs)
c lelli(i): logical flag: are crowns of class i ellipsoids?
c
c ------  CANOPY MEAN STRUCTURAL PARAMETERS (flag l_struc)
      double precision ulg(nclmax), uuu(nclmax), efflai, sntr, hmtree 
     & , vhekm, vhcilm, rmcrown, dbhmean, rmassm, slwm, vliit, cano
     & , tlty, tlaief, tlai, tbai   
      save ulg, uuu, efflai, sntr, hmtree, vhekm, 
     & vhcilm, rmcrown, dbhmean, rmassm, slwm, vliit, cano, tlty, 
     & tlaief, tlai, tbai   
c the above variables are described after the call to strmean() below
      double precision DIFN, optPAI, p_Pola
      save DIFN, optPAI, p_Pola
c DIFN: DIFfuse Non-interceptance, DIFN=1-i0
c optPAI: optical, effective LAI (from Miller's equation)
c p_Pola: Stenberg 2007, Eq.17 [Remote Sensing of Environment 109, 221–224]
c ---- gap probabilities in crown layer: unidirectional gaps
      double precision gaps_v(1), gaps_s(1), gaps_q(nknotm)
c     bidirectional gaps
      double precision psgvu(1)
      double precision bdgfu(nclmax) 
      double precision btr1uk(nclmax)
      double precision bdgfd(nclmax)
      double precision btr1dk(nclmax)
      save psgvu, bdgfu, btr1uk, bdgfd, btr1dk
c   gaps_v: gap probability, view direction
c   gaps_s: gap probability, Sun direction
c   gaps_q: gap probability, quadrature directions
c   psgvu: bidirectional (sun, view) gap probability, upper hemisphere
c      NOTE: gaps_v, gaps_a, psgvu are arrays of length 1 to make use of "sequence association" and pass it to 
c         a function excepting an array. This is cumbersome, but should be standards-compliant 
c   ===WARNING== the contents from the next variables were guessed from
c     comment lines without going through the code line by line!
c   bdgfu: Integral of bidirectional gap probability over a crown, upward direction
c   bdgfd: Integral of bidirectional gap probability over a crown, downward direction
c   btr1uk: probability of seeing sunlit trunk from above
c   btr1dk: probability of seeing sunlit trunk from below
c   The same bidirectional gap probabilities for quadrature
c   NOTE: compared with original frt, the dimensionality of the arrays is increased by one
c     and the dimensions are reordered. Transposed old variables can be obtained by fixing the 
c     last dimension (azimuth = phi). "Q" was added to the end of the variable names (for Quadrature).
      double precision psgvuQ(nknotm,nphim),
     & bdgfuQ(nclmax,nknotm,nphim), btr1ukQ(nclmax,nknotm,nphim),
     & bdgfdQ(nclmax,nknotm,nphim), btr1dkQ(nclmax,nknotm,nphim)
      save psgvuQ, bdgfuQ, bdgfdQ, btr1ukQ, btr1dkQ
c   psgvuQ: bidirectional (sun, view) gap probability, upper hemisphere
c   ===WARNING== the contents from the next variables were guessed from
c     comment lines without going through the code line by line!
c   bdgfuQ: Integral of bidirectional gap probability over a crown, upward direction
c   bdgfdQ: Integral of bidirectional gap probability over a crown, downward direction
c   btr1ukQ: probability of seeing sunlit trunk from above
c   btr1dkQ: probability of seeing sunlit trunk from below

c ------  OPTICAL PROPERTIES BY CLASS ------
      double precision rind, rnlfi, refl, tran
      double precision rlfcl(nclmax), tlfcl(nclmax), rnlf(nclmax),
     & rbrnc(nclmax), rtrnk(nclmax)
      save rlfcl, tlfcl, rnlf, rbrnc, rtrnk
      
c rind: leaf wax refractive index at current wavelength for current class 
c rnlfi: leaf wax reflractive index to be applied (after correction) at current wavelength for current class 
c refl: leaf reflectance at current wavelength for current class 
c tran: leaf transmittance at current wavelength for current class 
c rlfcl(): leaf diffuse (Lambertian) reflectances component for all tree classes at current wavelength
c tlfcl(): leaf transmittance for all tree classes at current wavelength
c rnlf(): leaf wax refractive index to be applied (after correction) for all tree classes at current wavelength
c rbrnc(): branch reflectance for all tree classes at current wavelength
c rtrnk(): trunk reflectance for all tree classes at current wavelength

c ------  CANOPY MEAN OPTICAL PROPERTIES, from optmean() ------
c     to be used in the  twostream submodel (recalculation flag l_opt .or. l_struc)
      double precision rleff, tleff, rneff, rty, pkhair,
     &  rrs(nclmax), ttt(nclmax), utot, rteff, tteff
      save rleff, tleff, rneff, rty, pkhair, rrs,
     & ttt, utot, rteff, tteff
c rleff: average reflectance of leaves+branches at current wavelength
c tleff: average tranmsmittance of leaves+branches at current wavelength
c rneff: average wax refractive index at current wavelength
c rty: average trunk reflectance  at current wavelength
c rrs(): average diffuse reflectace of a tree class, excl. trunks  at current wavelength
c ttt(): ??? at current wavelength
c utot: total effective PAI, leaves + branches + stems  at current wavelength
c rteff: average effective reflectance of canopy elements  at current wavelength
c tteff: average effective transmittance of canopy elements (incl. trunks) at current wavelength
c
c ------  FOREST FLOOR ------ 
c    filled from input data
      double precision rsogrou, rddgrou, rdsgrou, rsdgrou
      save rsogrou, rddgrou, rdsgrou, rsdgrou
c rsogrou: BRDF of ground vegetation (directional-directional reflectance) at current wavelength
c     NOTE: for flux computations, rsogrou is set to rsdgrou
c rddgrou: hemispherical-hemispherical forest floor reflectance (=albedo) at current wavelength
c rdsgrou: hemispherical-directional reflectance of ground vegetation at current wavelength

c ------  CUBATURE (QUADRATURE) ------ 
      integer nquad_t, nquad_p 
      double precision wght_q(nknotm), gq_w(nknotm),
     & xquad_t(nknotm), xquad_p(nphim)
      save nquad_t, nquad_p, wght_q, xquad_p, xquad_t
c  gq_w: Gauss quadrature weights, temporary variable used to calculate wght_q
c  nquad_t: number of knots over theta for canopy-level integration
c  nquad_p: number of knots over phi for canopy-level integration
c    Quadrature values and weights
c  xquad_t: quadrature nodes over theta (polar angle)
c  xquad_p: quadrature nodes over phi (azimuth angle)
c  wght_q: quadrature weights for canopy-level integration

c ------  RESULTS ------
c   ground + understory reflectances
c     most stored to SPCOUT at the end (depeding on ijob)
      double precision bi0u, bi0d, R1_c, R1_g, R1, T1, R, T
      double precision rhd_hi_c, rhd_hi_g, thd_hi

c bi0u: 1st order BRF (crowns only) at current wavelength
c bi0d: 1st order BTF (crowns only) at current wavelength
c R1_c: 1st order reflectance component (in view direction) from crowns
c R1_g: 1st order ground reflectance component
c R1: total first-order reflectance component in view direction
c T1: 1st order transmittance component 
c R: canopy hemispherical-directional reflectance factor
c T: canopy hemispherical-directional transmittance factor
c rhd_hi_c: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-direcitpnal), 
c     contribution by canopy, includes first-order diffuse-sky reflectance component
c rhd_hi_g: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-direcitpnal),
c     contribution by canopy, includes first-order diffuse-sky reflectance component
c thd_hi: higher-order component (excl. 1st-order direct) of canopy+soil hemispherical-directional transmittance, includes first-order diffuse-sky reflectance 
c NOTE:  the higher-order components (*_hi) include first-order reflectance of diffuse-sky radiation

c  integrated quantities for flux reflectance/transmittance
c     most stored to SPCOUT at the end (depeding on ijob)
      double precision f_up, f_up1_c, f_up1_g,
     &  f_uphi_c, f_uphi_g, f_down, f_down1

c  f_up: total flux reflectance
c  f_up1_c: 1st order direct flux reflectance component by crowns, f_up1_c is included in f_up
c  f_up1_g: 1st order direct flux reflectance component by ground, f_up1_g is included in f_up
c  f_uphi_c: higher order + 1st-order dif. sky flux reflectance component by crowns, f_up1_c is included in f_up
c  f_uphi_g: higher order + 1st-order dif. sky flux reflectance component by ground, f_up1_g is included in f_up
c  f_down: scattered (+diffuse) flux transmittance
c  f_down1: 1st order flux transmittance component, f_down1 is included in f_downs, total downward flux = f_downs + uncollided transmittance

c ------  MISC ------
      double precision c_fact
c     c_fact: correction factor for diffuse fluxes (to correct for the error in estimating 1st order scattering) at current wavelength
	  double precision CrC_i
c     CrC_i: temp variable for Crown Cover of a class

c     current wavelength 
      double precision rlambda, x1, x2
c 
c     counters etc.
      integer jcl, i, jwl, jqt, jqa, jth
c     control flags
      logical l_opt, l_struc, l_dirv, l_dirs, l_grnd, l_refl
c      A flag for whether we are called for flux computations
      logical l_flux
c     how many zenith and azimuth angles are used 
      integer n_theta, n_phi
      save n_theta, n_phi

c --------------------- start the work

c     this is done also in frt main program, but in case comprt is called from somewhere else
c     /pidr/ is not used in comprt, but in some subroutines called by it
c     at some point, /pidr/ should be made redundant
      pi = acos(-1.d0)
      dr = pi/180.d0
      
      l_opt =    ( IC(4) .gt. 0 )
      l_struc =  ( IC(5) .gt. 0 )
      l_dirv =   ( IC(6) .gt. 0 )
      l_dirs =   ( IC(7) .gt. 0 )
      l_grnd =   ( IC(8) .gt. 0 )
      l_refl =   ( IC(9) .gt. 0 )
      l_flux = .false.
      if ( IC(12+nclmax) .eq. -1 ) l_flux = .true.
      
      if ( l_struc ) then
        write(*,*)'comprt: computing structural parameters'
        nzs = IC(3)
c       Compute knots and weights of Gauss-Legendre quadrature over the hemisphere
c       nquad_t and nquad_p are the number of knots
c          and xquad_t() and xquad_p() for the knots themselves
c       Azimuth: the model is symmetric wrt principal plane, divide nodes equally in [0,pi]      
c          if nquad_p == 1, this leads to the cross plane      
        nquad_t = IC(13+nclmax)
        nquad_p = IC(14+nclmax)
        do i = 1, nquad_p
          xquad_p(i) = pi*(1.0/nquad_p/2.0 + (i-1.0)/nquad_p)
        enddo
c       The quadrature weights combine both theta and phi weights so that cosine
c       -weighed integral of f(x) over a hemisphere is 
c         F = SUM_i[ w_i * SUM_j( f(theta_i,phi_j) ) ]
        call gauleg(0.d0, pi/2.d0, xquad_t, gq_w, nquad_t)
        do jth = 1, nquad_t
          wght_q(jth)=2*gq_w(jth)*cos(xquad_t(jth))*sin(xquad_t(jth))
     &        / nquad_p
        enddo
        if ( .not. l_volint ) then
          write(*,*)'comprt: filling internal cubature'
c         Fill the knots and weights of cubature on a sphere (for integrating
c         gap probability over a crown ellipsoid)
          call cubell9(netst, xetst, yetst, zetst, aetst)
c         Knots and weights of cubature on a circle (for integrating
c         gap probability over a cylindrical crown)
          call cubcirc(nctst, xctst, yctst, actst)
c         Calculate the points of the quadrature for the vertical direction
          x1      = -1.d0
          x2      = 1.d0
          call gauleg(x1, x2, zctst, acztst, 2*nzs)
          l_volint = .true.
        endif
c       for clarity, copy the variables from the input vector XC
c       into better-named variables at first call of the function
        ncl = IC(1)
        do jcl = 1, nclmax
          lelli(jcl) = ( IC(jcl+11) .gt. 0 ) ! are crowns ellipsoids?
c          lelli = ( IC(i+11) > 0, i = 1, ncl )
          stdns(jcl)  = XC(jcl, 1) ! stand density [m^-2]
          htr(jcl)    = XC(jcl, 2) ! tree height [m]
          hc1(jcl)    = XC(jcl, 3) ! crown length, ell / con [m]
          hc2(jcl)    = XC(jcl, 4) ! crown length, cylinder [m]
          if (hc1(jcl) .gt. htr(jcl)) then
            hc1(jcl) = htr(jcl)
            print *, 'FUNC.F:'
            print *, 'Warning: tree height is less than crown length'
            print *, jcl, ' -> crown length is adjusted.'
          endif
          rcr(jcl)    = XC(jcl, 5) ! crown radius [m]
          dbh(jcl)    = XC(jcl, 6)*1.d-2 ! trunk diameter [m]
          rmass(jcl)  = XC(jcl, 7) ! dry leaf wght [kg/tree]
          slwcl(jcl)  = XC(jcl, 8) ! SLW [gm^-2]
          rlai(jcl)   = stdns(jcl)*rmass(jcl)*1.d3/slwcl(jcl) ! LAI
          rbai(jcl)   = XC(jcl, 9)*rlai(jcl) ! BAI
          clmpst(jcl) = XC(jcl, 10) ! tree distr. parameter
c                                      =cB in eq. (11) in Nilson 1999
c         calculate glmp: Fisher's grouping index
          call iterats(clmpst(jcl), glmp(jcl))
c         correct for the unrealistic case of canopy cover  > crown cover by adjusting Fisher's index FGI
c         according to Nilson and Kuusk (2004, Agricultural and Forest Meteorology 124, 157–169)
c         if no overlapping crowns appears (canopy cover = crown cover), FGI = 1-CrownCover.
c         This sets the practical lower limit on the FGI regularity: even lower FGI would indicate even more
c         regular distribution, but because a limit has been reached, this would not decrease canopy transmittance.
          CrC_i = stdns(jcl)*rcr(jcl)**2 
          if ( glmp(jcl) .lt. (1.0-CrC_i) ) then
            glmp(jcl) = 1.0-CrC_i
            write(*,'(a,i2,a,F5.3)')
     &        "corrected FGI of tree class",jcl,":",glmp(jcl)
		  endif
          clmpsh(jcl) = XC(jcl, 11) ! shoot shading coef (ssc)
          crncl(jcl) = XC(jcl, 12) ! wax refration index ratio
          shl(jcl)   = XC(jcl, 13)
c       NOTE: in original frt, crncl and shl were at position 25 and 26, respectively
          if (shl(jcl) .le. 0.d0) shl(jcl) = .01d0
c WWW: warn when modifying sth! 
        enddo
      endif
      if (l_refl ) then
        if ( l_dirs ) then 
          thetv = XC(1,21)
          phiv = XC(1,22)
          if ( l_flux ) then
            n_theta =  nquad_t
            n_phi = nquad_p
          else
            n_theta = 1
            n_phi = 1
          endif
        endif
      else
!       reflectance factors not computed, make sure loops over angles will not be run
        n_theta = 0
        n_phi = 0
        l_flux = .false.
      endif
      if ( l_dirs ) then
        thets = XC(1,23)
      endif 
      pkhair = XC(1,24)
c
c  CALCULATE MEAN PARAMETERS FOR THE STAND (averaged over tree classes)
c
        if ( l_struc ) then
          call strmean
     &    (lelli, ncl, stdns, htr, hc1, hc2, rcr, dbh,
     &     rmass, slwcl, rlai, rbai, clmpst, clmpsh, glmp, ulg, uuu,
     & efflai, sntr, hmtree, vhekm, vhcilm, rmcrown, dbhmean, rmassm,
     &     slwm, vliit, cano, tlty, tlaief, tlai, tbai)
c         strmean outputs:
c         ulg: uuu / 2 ???
c         uuu: effective plant area density (leaves corrected for clumping + branches)
c         efflai: LAI of a random canopy that causes extinction similar to that of the stand at thseff (=40°)
c         sntr: total number of trees per m^2
c         hmtree: mean tree height
c         vhekm: mean crown length (ellipsoid/cone)
c         vhcilm: mean crown length (cylinder)
c         rmcrown: mean crown radius
c         dbhmean: mean dbh
c         rmassm: mean leaf mass per m^2
c         slwm: mean SLW
c         vliit: crown cover
c         cano: canopy cover
c         tlty: total unshaded trunk area per m^2 (below crown, calculated as pi*dbh*l0)
c         tlaief: total effective LAI (corrected for shoot-level clumping)
c         tlai: total lai (leaves only)
c         tbai: total BAI (branch area index) so that PAI = tlai+tbai
c
c         initialize structural variables, unidirectional gap probabilities for quadrature
c
          call hetk8sA
     &     (lelli, ncl, shl, nzs, nquad_t, xquad_t, 
     &      stdns, htr, hc1, hc2, rcr, dbh, glmp, ulg, uuu, gaps_q)
c          calculate also canopy recollision probability and effective (optical) LAI
c           using Stenberg (2007) and Miller (1967), respectively
          DIFN = 0
c         DIFfuse Non-interceptance DIFN=1-i0
          optPAI = 0
          do jqt = 1,nquad_t
            DIFN = DIFN + wght_q(jqt)*gaps_q(jqt)*nquad_p
c           multiplication by nquad_p simulates the different azimuth angles in the quadrature,
c               each giving the same result
            optPAI = optPAI - wght_q(jqt)*log(gaps_q(jqt))*nquad_p
          enddo
c         note: p is calculated using PAI, not just LAI
          p_Pola = 1-(1-DIFN)/(tlai+tbai)
c           Stenberg 2007, Eq.17 [Remote Sensing of Environment 109, 221–224]
          XOUT(1) = sntr
          XOUT(2) = hmtree
          XOUT(3) = vhekm
          XOUT(4) = vhcilm
          XOUT(5) = rmcrown
          XOUT(6) = dbhmean
          XOUT(7) = rmassm
          XOUT(8) = slwm
          XOUT(9) = tbai
          XOUT(10) = tlaief
          XOUT(11) = tlai
          XOUT(12) = vliit
          XOUT(13) = cano
          XOUT(14) = DIFN
          XOUT(15) = optPAI
          XOUT(16) = p_Pola
          XOUT(20) = efflai
        endif
c       precompute bidirectional gap probabilities for view and sun angles and quadrature
c       this avoids calling hetk8sB separately for each direction
        if (l_struc .or. l_dirs .or. l_dirv ) then
          if ( l_flux ) then
c           flux calculations, quadrature used
c           loop over azimuth, hetk8sB contains internal loop over zenith angle
c           if no flux computations, hetk8sB is called from within wavelength loop
            do jqa = 1, nquad_p
              call hetk8sB
     &    (lelli, ncl, shl, nzs, nquad_t, xquad_t, thets, xquad_p(jqa),
     &      stdns, htr, hc1, hc2, rcr, dbh, glmp, ulg, uuu,
     &      psgvuQ(1,jqa), bdgfuQ(1,1,jqa), bdgfdQ(1,1,jqa),
     &      btr1ukQ(1,1,jqa), btr1dkQ(1,1,jqa) )
            enddo
          endif
        endif
        if ( l_dirs .or. l_struc ) then
c          calculate gaps_s
          call hetk8sA
     &     (lelli, ncl, shl, nzs, 1, [thets], 
     &      stdns, htr, hc1, hc2, rcr, dbh, glmp, ulg, uuu, gaps_s) 
        endif
c  
c     ------------------------------------------------------------------
c                                                  LOOP OVER WAVELENGTHS
      write (*,*) "comprt starting loop over wavelengths..."
      do jwl = IC(10), IC(11)
        rlambda = SPCIN(jwl,1)
        write(*,'(a16,i5,f8.1)')'wavelength #',jwl,rlambda
c         rlambda -- wavelength (nm)
c         jwl -- index in all possible wavelengths (starting at 400 nm with 1 nm step)

        sqratio = SPCIN( jwl, 6 )
        c_fact = SPCIN( jwl, 8+5*nclmax )
c
c       LEAF, BRANCH, TRUNK REFLECTANCE AND TRANSMITTANCE
c
        if ( l_opt) then
c          write(*,*)'DBG: comprt recalculating optical properties'
          rind = SPCIN(jwl,7) ! leaf wax refractive index
          do jcl = 1, ncl
            rnlfi      = rind*crncl(jcl)
            if (rnlfi .lt. 1.d0) rnlfi = 1.d0
c             setting crncl to, e.g., 0 can thus be used to ignore specular reflectance
c             (if rnlfi=1, refractive indices are equal and specular reflectance is 0)

            refl = SPCIN(jwl,7+jcl)
            tran = SPCIN(jwl,7+nclmax+jcl)         
c           separate diffuse and specular reflectance components
c               specular component calculated for normal incidence
            rlfcl(jcl) = refl - ((1.d0 - rnlfi)/(1.d0 + rnlfi))**2
            tlfcl(jcl) = tran
            rnlf(jcl)  = rnlfi
            rbrnc(jcl) = SPCIN(jwl,7+3*nclmax+jcl)
            rtrnk(jcl) = SPCIN(jwl,7+4*nclmax+jcl)
          enddo
        endif ! if ( l_opt) then

        if ( l_opt .or. l_struc ) then
c         The mean values of optical parameters for computing diffuse 
c           reflectance with the  2-stream submodel
c          write(*,*)'DBG: comprt recalculating 2-stream input'
          call optmean
     &    ( ncl, stdns, htr, hc1, hc2, dbh, rlai, rbai, clmpsh, tlty,
     &      tlaief, tbai, rlfcl, tlfcl, rnlf, rbrnc, rtrnk,
     &      rleff, tleff, rneff, rty, rrs, ttt, utot, rteff, tteff )
c          rleff: average reflectance of leaves+branches [not used?]
c          tleff: average tranmsmittance of leaves+branches [not used?]
c          rneff: average corrected wax refractive index
c          rty: average trunk reflectance
c          rrs: average diffuse reflectance of a tree class, excl. trunks
c          ttt: ???
c          utot: total effective LAI, leaves + branches + stems (not used?)
c          rteff: average effective reflectance of canopy elements (incl. trunks)
c          tteff: average effective transmittance of canopy elements (incl. trunks)
        endif
c       Initialize for flux calculation if needed
        if ( l_flux ) then
          write(*,*) 'comprt initializing flux computations'
          f_up = 0 
          f_up1_c = 0 
          f_up1_g = 0 
          f_uphi_c = 0
          f_uphi_g = 0
          f_down = 0 
          f_down1 = 0 
        endif
c ===== outer loop over view zenith ========================================
        do jqt = 1, n_theta
          if ( l_flux ) then
            thetv = xquad_t(jqt)
            l_dirv = .true. ! set to .true. in the beginning of each loop
          endif
c ====== inner loop over view azimuth  =====================================
          do jqa = 1, n_phi
            if ( l_flux ) then
              phiv = xquad_p(jqa)
              l_dirv = .true. ! set to .true. in the beginning of each loop
            endif
c
c           CALCULATE OR LOAD GAP PROBABILITIES IN CROWN LAYER
c
            if ( l_dirv .or. l_dirs .or. l_struc ) then
              if ( l_flux ) then
c               load gaps_v and bidirectional variables
                gaps_v(1) = gaps_q( jqt )
                psgvu(1) = psgvuQ( jqt, jqa )
                do jcl = 1,ncl
                  bdgfu(jcl) = bdgfuQ( jcl, jqt, jqa )
                  btr1uk(jcl) = btr1ukQ( jcl, jqt, jqa )
                  bdgfd(jcl) = bdgfdQ( jcl, jqt, jqa )
                  btr1dk(jcl) = btr1dkQ( jcl, jqt, jqa )
                enddo
              else
c             compute bidirectional gap probabilities
                call hetk8sB
     &          (lelli, ncl, shl, nzs, 1, [thetv], thets, phiv,
     &          stdns, htr, hc1, hc2, rcr, dbh, glmp, ulg, uuu,
     &          psgvu, bdgfu, bdgfd, btr1uk, btr1dk )
                if ( l_dirv ) then
c                 calculate gaps_v
                  call hetk8sA
     &            (lelli, ncl, shl, nzs, 1, [thetv], stdns, htr,
     &             hc1, hc2, rcr, dbh, glmp, ulg, uuu, gaps_v)
                endif
              endif
c             all structural variables computed
              l_struc = .false.
            endif
c           LOAD GROUND REFLECTANCES
c             rsogrou: BRDF of ground vegetation (directional-directional reflectance)
c               NOTE: for flux computations, rsogrou is not used, rsdgrou is used instead
c             rsdgrou: directional-hemispherical forest floor reflectance
c             rddgrou: hemispherical-hemispherical forest floor reflectance (=albedo)
c             rdsgrou: hemispherical-directional reflectance of ground vegetation
            if ( l_grnd ) then
              rddgrou = SPCIN( jwl , 2 )
              rsdgrou = SPCIN( jwl , 3 )
              rdsgrou = SPCIN( jwl , 4 )
              if ( l_flux ) then
                rsogrou = rdsgrou
              else 
                rsogrou = SPCIN( jwl , 5 )
              endif
            else
c             directional components are reloaded when required by the specific flags
c             WWW rsogrou and rdsgrou are reloaded in the flux loop every time, although they do not change.
              if ( l_dirs .or. l_dirv ) then
                if ( l_flux ) then
                  rsogrou = SPCIN( jwl , 4 )
                else 
                  rsogrou = SPCIN( jwl , 5 )
                endif
              endif
              if ( l_dirv ) then
                rdsgrou = SPCIN( jwl , 4 )
              endif
              if ( l_dirs ) then
                rsdgrou = SPCIN( jwl , 3 )
              endif
            endif
c           nothing specific to solar and view directions left in the loops
            l_dirs = .false.
            l_dirv = .false.

c           ACTUAL CALCULATION OF REFLECTANCE OF TRANSMITTANCE
c           === het8o and twostr inputs ===
c           ncl: no. of tree classes
c           phiv: relative view azimuth, rad
c           stdns: stand density
c           uuu: effective plant area density (leaves corrected for clumping + branches)[strmean] vector(nclmax)
c           efflai: LAI of a random canopy that causes
c             extinction similar to  that of the stand at thseff (= 40°) [strmean]
c           sqratio: direct to total irradiance at TOC
c           rteff: average effective reflectance of canopy elements [optmean]
c           tteff: ??? average leaf diffuse reflectance (incorrectly calculated)? [optmean]
c           tlty: total unshaded trunk area per m^2 (below crown, calculated as pi*dbh*l0) [strmean]
c           rtrnk : branch reflectance at current wl, vector(nclmax) [SPCIN]
c           c_fact: correction factor for diffuse fluxes [SPCIN]
c           bdgfu: Integral of bidirectional gap probability over a crown, upward direction [hetk8sB]  vector (nclmax)
c           bdgfd: Integral of bidirectional gap probability over a crown, downward direction [hetk8sB]  vector (nclmax)
c           btr1uk: probability of seeing sunlit trunk from above [hetk8sB] (nclmax)
c           btr1dk: probability of seeing sunlit trunk from below [hetk8sB] vector (nclmax)
c           pkhair: leaf hair optical index (vaguely defined, set to 0.1), input parameter XC(1,24)
c                        pkhair affects specular gamma-function in gmfres
c           rnlf: wax refractive index, calculated in comprt.f
c           rrs: average diffuse reflectace of a tree class, excl. trunks [optmean]
c           ttt: mean element spectral transmittance for each tree class [optmean]
c           rsogrou: BRDF of ground vegetation (directional-directional reflectance) [SPCIN]
c           rdsgrou: hemispherical-directional reflectance of ground vegetation [SPCIN]
c           rsdgrou: directional-hemispherical reflectance of ground vegetation [SPCIN]
c           rddgrou: hemispherical-hemispherical ground reflectance (=albedo) [SPCIN]
c
c           calculate single scattering from tree crowns
            call hetk8o     
     &        (ncl, thets, thetv, phiv, stdns, uuu,
     &          bdgfu, bdgfd, btr1uk, btr1dk, pkhair, rnlf, 
     &          rtrnk, rrs, ttt, bi0u, bi0d )
c          == hetk8o outputs ==
c             bi0u: single scattering reflectance component (tree crowns)
c             bi0d: single scattering transmittance component (tree crowns)
 
            R1_c = sqratio*bi0u !   1st order reflectance component (in view direction) from crowns
            R1_g = sqratio*rsogrou*psgvu(1) !  1st order ground reflectance component
            R1 = R1_c + R1_g ! total first-order reflectance component in view direction
            T1 = sqratio*bi0d !  1st order transmittance component 
     
            call twostr
     &       ( thets, thetv, sqratio, efflai, tlty,
     &        c_fact, rteff, tteff, pkhair, rsogrou, rdsgrou,
     &        rsdgrou, rddgrou, gaps_s(1), gaps_v(1),
     &        rhd_hi_c, rhd_hi_g, thd_hi )
c          WWW rsogrou not used in twostr
c          == twostr outputs ===
*           rhd_hi_c: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-direcitpnal), 
*               contribution by canopy, includes all-order diffuse-sky reflectance component
*           rhd_hi_g: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-direcitpnal),
*               contribution by ground (forest floor), includes all-order diffuse-sky reflectance component
!           thd_hi: higher-order component (excl. 1st-order direct) of canopy+soil hemispherical-directional transmittance 
*               includes all-order diffuse-sky reflectance component
*           the higher-order components (*_hi) include all orders of reflectance of diffuse-sky radiation
*           NOTE:  the higher-order components (*_hi) include
*                *  first-order reflectance of diffuse-sky radiation weighted by diffuse sky flux contribution
*                *  higher order reflectance of direct sun and diffuse sky fluxes, 
*                      weighted by direct and diffuse sky TOC flux contributions, respectively.

            R = R1 + rhd_hi_c + rhd_hi_g ! total forest reflectance in view direction
            T = T1 + thd_hi

            if ( l_flux ) then
c             add contribution to fluxes
              f_up = f_up + R*wght_q(jqt)
              f_up1_c = f_up1_c + R1_c*wght_q(jqt)
              f_up1_g = f_up1_g + R1_g*wght_q(jqt)
              f_uphi_c = f_uphi_c + rhd_hi_c*wght_q(jqt)
              f_uphi_g = f_uphi_g + rhd_hi_g*wght_q(jqt)
              f_down = f_down + T*wght_q(jqt)
              f_down1 = f_down1 + T1*wght_q(jqt)
            endif
c           Loops over geometry end here
          enddo
c            do jqa = 1,nquad_p
        enddo
c          do jqt = 1, n_theta
c       store output
        if (l_flux) then
          SPCOUT(jwl, 1) = f_up
          SPCOUT(jwl, 2) = f_down
          SPCOUT(jwl, 3) = f_up1_c
          SPCOUT(jwl, 4) = f_up1_g
          SPCOUT(jwl, 5) = f_uphi_c
          SPCOUT(jwl, 6) = f_uphi_g
          SPCOUT(jwl, 7) = f_down1
c         when computing fluxes, the gap probabilities involving view direction are rather pointless
          XOUT(17) = DIFN
          XOUT(19) = DIFN*gaps_s(1)
        else
          SPCOUT(jwl, 1) = R 
          SPCOUT(jwl, 2) = T 
          SPCOUT(jwl, 3) = R1_c
          SPCOUT(jwl, 4) = R1_g
          SPCOUT(jwl, 5) = rhd_hi_c
          SPCOUT(jwl, 6) = rhd_hi_g
          SPCOUT(jwl, 7) = T1
          XOUT(17) = gaps_v(1)
          XOUT(19) = psgvu(1)
        endif
        SPCOUT(jwl, 8) = rteff
        SPCOUT(jwl, 9) = tteff
c       store final missing elements in XOUT
        XOUT(18) = gaps_s(1)
c       for next wavelength value: set flags for recalculations of params
        l_grnd = .true.
        l_opt = .true.
      enddo 
c     =========== END OF LOOP OVER WAVELENGTHS ==========================

      return
      end
      
