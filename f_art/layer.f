!   part of the fortran-only version of FRT
	  subroutine layer
     & (thets, thetv, ggs, ggv, ggj, ul, rrl, ttl,
     & rdd, tdd, rdo, tdo, rsd, tsd, rso)
c
c  Diffuse reflection and transmission of a layer from the 2-stream model
c  for elliptical leaf orientation. A. Kuusk  13.01.2000
*  Calculates the scattering operators of a leaf layer (Table 1; Kuusk 2001)
*  Code is similar to that given in Kuusk 1995, A computer-efficient...
*     Comp & Geosci 22, 149-163, Appendix B
*  and Kuusk 1995, A fast invertible..., RSE 51, 342-350
*
* IN:
c  thets - the Sun zenith, thetv - view angle
c  ggs, ggv - G-function at thets and thetv
c  ggj  - the J-integral divided by 2, eq. (6)
c  ul - LAI
c  rrl, ttl - leaf reflection and transmission
* OUT:
c  Notation: [r=reflectance,t=transmittance][s=sun,d=diffuse][o=observer,d=diffuse]
c  rso - higher-order reflectance of the layer in view direction
c        for direct incoming radiation
c  Other outputs include most probably all scattering orders, but this is contradictory in (Kuusk, 2001)
*    Kuusk (2001) Eq. (9) for rho_d^plants includes SQr_so, although it should not include by definition
*    first-order scattering from the direct beam. Exclusion of first-order scattering from r_so is not
*    explicitly stated anywhere, but is evident from Eq. (9).

      implicit double precision (a-h, o-z)
      include 'frtpar.h'
c
*     eint(xy) = (1.d0 - exp(-xy))/xy
c
      rtp   = (rrl + ttl)/2.d0
      rtm   = (rrl - ttl)/2.d0
      bf    = rtm*ggj
c
*  Eq. numbers in Kuusk 1996, A computer-efficient, Comp & Geosci 22, 149-163, App. B
      vks   = ggs/cos(thets)        ! k_1 from B30, excl u_L
      vkv   = ggv/cos(thetv)        ! k_2 from B30, excl u_L
      vsig  = rtp + bf                       ! B23, excl u_L
      vatt  = 1.d0 - rtp + bf                ! B22, excl u_L
      vssf  = vks*rtp - bf                   ! B24, excl u_L
      vssb  = vks*rtp + bf                   ! B25
      vmm   = sqrt(vatt**2 - vsig**2)        ! B21
      vvv   = vkv*rtp + bf
      vuu   = vkv*rtp - bf
      vh1   = (vatt + vmm)/vsig              ! B20
      vh2   = 1.d0/vh1
      vcc   = (vssf*vsig - vssb*(vks - vatt))/(vmm**2 - vks**2) ! B16
      vdd   = (vssb*vsig + vssf*(vks + vatt))/(vmm**2 - vks**2)
c
      exmu1 = exp(vmm*ul) ! ???
      exmu2 = 1.d0/exmu1
      exku1 = exp(-vks*ul)
      delt  = vh1*exmu1 - vh2*exmu2
c
      rdd   = (exmu1 - exmu2)/delt
      tdd   = (vh1 - vh2)/delt
c
      delta = vh2*vcc*exku1 - vdd*exmu1
      aa1   = delta/delt
      deltb = vdd*exmu2 - vh1*vcc*exku1
      bb1   = deltb/delt
      tsd   = vh1*aa1*exmu2 + vh2*bb1*exmu1 + vdd*exku1
      rsd   = aa1 + bb1 + vcc
      tdo   = ul*((vuu*vh1 + vvv)*eint((vkv - vmm)*ul) -
     &        (vuu*vh2 + vvv)*eint((vkv + vmm)*ul))/delt
      rdo   = ul*(exmu1*(vvv*vh1 + vuu)*eint((vkv + vmm)*ul) -
     &        exmu2*(vvv*vh2 + vuu)*eint((vkv - vmm)*ul))/delt
c
      rso   = ul*((vvv*vdd + vuu*vcc)*eint((vkv + vks)*ul) +
     &        (vvv*vh1 + vuu)*aa1*eint((vkv + vmm)*ul) +
     &        (vvv*vh2 + vuu)*bb1*eint((vkv - vmm)*ul))
c

      return
      end
*
**********************************************************************
*
