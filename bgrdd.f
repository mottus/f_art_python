!   part of the fortran-only version of FRT
      subroutine bgrdd
     & (thets, thetv, skyl, efflai, tlty,
     & rteff, tteff, rddgrou, rsdgrou, bddif)
c   diffuse fluxes down
c                                  A. Kuusk  14.10.2000, 12.09.2002
c   tlty - LAI of trunks
c
      implicit double precision (a-h, o-z)
      include 'frtpar.h'

      data gsf/.5d0/
c
*     eint(xy) = (1.d0 - exp(-xy))/xy
c
c                                      Difuussed vood 2-voo lähenduses
c
      ul    = efflai + tlty
      gcyl  = 2.d0*sin(thets)/pi
      ggef  = (efflai*gsf + tlty*gcyl)/ul
      ggs   = ggef
      ggv   = ggef
      ggj   = efflai/6.d0/ul
c
      rtp  = (rteff + tteff)/2.d0
      rtm  = (rteff - tteff)/2.d0
      bf   = rtm*ggj
c
*     cthetv = cos(thetv)
c                                   ! v??? - Verhoefi tähistused
      vks   = ggs/cos(thets)
      vkv   = ggv/cos(thetv)
      vsig  = rtp  + bf
      vatt  = 1.d0 - rtp + bf
      vssf  = vks*rtp - bf
      vssb  = vks*rtp + bf
      vmm   = sqrt(vatt**2 - vsig**2)
      vvv   = vkv*rtp + bf
      vuu   = vkv*rtp - bf
      vh1   = (vatt + vmm)/vsig
      vh2   = 1.d0/vh1
      vcc   = (vssf*vsig - vssb*(vks - vatt))/(vmm**2 - vks**2)
      vdd   = (vssb*vsig + vssf*(vks + vatt))/(vmm**2 - vks**2)
c
      exmu1 = exp(vmm*ul)
      exmu2 = 1.d0/exmu1
      exku1 = exp(-vks*ul)
c
      delt  = exmu2*(vh2 - rddgrou) - (vh1 - rddgrou)*exmu1
      delta = (1.d0 - skyl)*exku1*(vdd*rddgrou + rsdgrou - vcc)*vh2 -
     &        (1.d0 - rddgrou*vh2)*(skyl - vdd*(1.d0 - skyl))*exmu1
      deltb = exmu2*(1.d0 - vh1*rddgrou)*(skyl - vdd*(1.d0 - skyl)) -
     &        vh1*(1.d0 - skyl)*exku1*(vdd*rddgrou + rsdgrou - vcc)
      vaa   = delta/delt
      vbb   = deltb/delt
c
      bddif = vaa*(vuu*vh1 + vvv)*ul*exmu2*eint(ul*(vkv - vmm)) +
     &        vbb*(vuu*vh2 + vvv)*ul*eint(ul*(vkv + vmm))*exmu1 +
     &        (vuu*vdd + vvv*vcc)*ul*(1.d0 - skyl)*exku1*
     &        eint(ul*(vkv - vks))
c
      return
      end
*
************************************************************************
*
