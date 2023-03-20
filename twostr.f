      subroutine twostr
     &  (thets, thetv, sqratio, efflai, tlty,
     &   cf, rteff, tteff, p_khair, rsogrou, rdsgrou,
     &   rsdgrou, rddgrou, t0, tv,
     &   Rhd_hi_c, Rhd_hi_g, Thd_hi )

!  Two-stream calculations for calculating diffuse fluxes
!
!  IN:
!  thets: solar zenith, rad
!  thetv: view zenith, rad
!  efflai: strmean
!  sqratio: SPCIN
!  rteff: optmean
!  tteff: optmean
!  tlty: strmean
!  rtrnk: SPCIN
!  cf: SPCIN
!  p_khair: leaf hair optical index (vaguely defined, set to unity), input parameter XC(1,24) NOT USED?
!  rsogrou: SPCIN -- NOT USED?
!  rdsgrou: SPCIN
!  rsdgrou: SPCIN
!  rddgrou: SPCIN
!  t0: gaps_s(1) from hetk8s: transmittance in solar direction
c  tv: gaps_s(1) from hetk8s: transmittance in view direction
!  OUT:
*  Rhd_hi_c: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-direcitpnal), 
*          contribution by canopy, includes all-order diffuse-sky reflectance component
*  Rhd_hi_g: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-direcitpnal),
*          contribution by ground (forest floor), includes all-order diffuse-sky reflectance component
!  Thd_hi: higher-order component (excl. 1st-order direct) of canopy+soil hemispherical-directional transmittance 
*          includes all-order diffuse-sky reflectance component
*    the higher-order components (*_hi) include all orders of reflectance of diffuse-sky radiation

      implicit none
      include 'frtpar.h'

      double precision thets, thetv, sqratio, efflai, tlty
      double precision cf, rteff, tteff, p_khair
      double precision rsogrou, rdsgrou, rsdgrou, rddgrou, t0, tv
      double precision Rhd_hi_c, Rhd_hi_g, Thd_hi
*    Notation: [reflectance,transmittance][hemispherical,directional,isotropic][directional,hemispherical]_[higherorder,total][Blacksurface]_[canopy,ground]
*    tdi_t: total directional-isotropic layer scattering factor (tsd in Kuusk 2001) FOR BLACK SOIL
*    rid_t: total isotropic-directional layer scattering factor (rsd Kuusk 2001) FOR BLACK SOIL
*    rii_t: total isotropic-isotropic layer scattering factor (rdd Kuusk 2001) FOR BLACK SOIL
*    tid_t: total isotropic-directional layer transmittance factor (tsd Kuusk 2001) FOR BLACK SOIL
*    tii_t: total isotropic-isotropic layer transmittance factor (tdd Kuusk 2001) FOR BLACK SOIL

      double precision tdi_t, rid_t, rii_t, tii_t, tid_t,
     &  rdi_t, rdd_hi, tdd_hi
     
      double precision gsf, utot, gcyl, ggef, ggj, fac, skyl

*     spherical G
      data gsf/.5d0/ 
  
      skyl  = 1.d0 - sqratio 
          
c     Compute layer scattering factors (black ground) 
      utot = efflai + tlty ! lai + trunk area index
      gcyl = 2.d0*sin(thets)/pi
      ggef = (efflai*gsf + tlty*gcyl)/utot
      ggj  = efflai/6.d0/utot
c     ggj  - the J-integral in the SAIL model * 0.5
      call layer
     & (thets, thetv, ggef, ggef, ggj, utot, rteff, tteff,
     & rii_t, tii_t, rid_t, tid_t, rdi_t, tdi_t, rdd_hi)
c     NOTE: layer does not output tsot, sun->observer diffuse transmittance.
c     NOTE: rdi_t not used below
      
c     Compute R,T which includes the effect of reflective ground 
c     correction factor is only applied to directional->X quantities
      fac  = 1.d0/(1.d0 - rii_t*rddgrou) 

      Rhd_hi_c   = sqratio*cf*rdd_hi + skyl*rid_t
     &       + (sqratio*t0*rsdgrou + sqratio*cf*tdi_t*rddgrou
     &       + skyl*tii_t*rddgrou)*tid_t*fac        ! Kuusk 2001,  (9)
c    
c      Rhd_hi_g = (sqratio*tdi_t + skyl*tii_t)*rdsgrou*tv*fac
c            ! Kuusk 2001, (10) with added fac -- probably en error in the original eqn
c     Modified in 2022 Rhd_hi_g to include the contribution by direct 
c       transmittance + multiple scattering 
       Rhd_hi_g=(sqratio*cf*tdi_t + skyl*tii_t + 
     &  sqratio*t0*rsdgrou*rii_t) *rdsgrou*tv*fac

c     the index _t added for clarity to indicate that the factor are for all scattering orders (total)
c     Compute higher-order down in the direction of observation,
c     AKA higher-order directional transmittance, Thd_hi
c     bgrff can output Thd_hi for a reflective ground, but its output has not been corrected 
c       with cf  and the subroutine is difficult to modify.
c       Hence, use bgrdd to compute tdd_hi (a.k.a. tsot in Kuusk 2001), the layer scattering
c       coefficient for higher-order scattering, and add the effect of reflecting surface later.
c       To compute tdd_hi, set skyl and ground reflectance to zero
      call bgrdd
     &     (thets, thetv, 0, efflai, tlty,
     &     rteff, tteff, 0, 0, tdd_hi)
c       compose Thd_hi from the layer scattering factors
        Thd_hi = sqratio*cf*tdd_hi + skyl*tid_t +
     &    ( sqratio*(t0*rsdgrou+tdi_t*cf*rddgrou) + 
     &    skyl*tii_t*rddgrou ) *rid_t / (1-rii_t*rddgrou)

c       If there is no need to correct, bgrdd can be used directly
c        call bgrdd
c     &    (thets, thetv, skyl, efflai, tlty,
c     &     rteff, tteff, rddgrou, rsdgrou, Thd_hi)

      return
      end
