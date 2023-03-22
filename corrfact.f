!    part of the frt distribution
!    Not used by the python version of FRT, python bindings for frt_wrapper.py
	  subroutine corrfact( XC, N_XC, IC, N_IC, SPCIN, N_SPCIN, 
     & N_XOUT, N_SPCOUT, cf )
!f2py intent(in) XC, N_XC, IC, N_IC, SPCIN, N_SPCIN, N_XOUT, N_SPCOUT
!f2py intent(out) cf
!  depends on comprt.o and all dependencies therein
!  compile with sth like
! f2py.exe -c --compiler=mingw32 -m corrfact comprt.o bck3.o bgrdd.o enel3.o hetk8o.o hetk8s.o layer.o optmean.o rmsub.o spooi.o strmean.o twostr.o corrfact.f
!     Compute the correctcompion factor for matching frt first- and higher-order scattering
!     See documentation and Mõttus, M., Stenberg, P. & Rautiainen, M. (2007). Photon recollision 
!       probability in heterogeneous forest canopies: compatibility with a hybrid GO model. 
!       Journal of Geophysical Research – Atmospheres 112, D03104, doi: 10.1029/2006JD007445.
      implicit none
      include 'frtpar.h'

      integer N_XC, N_IC, N_SPCIN, N_XOUT, N_SPCOUT
      double precision XC(nclmax,N_XC), SPCIN(nspchnl,N_SPCIN)
      double precision cf(nspchnl)
      integer IC(N_IC)
      
      double precision XC2(nclmax,N_XC), SPCIN2(nspchnl,N_SPCIN)
      double precision XOUT(N_XOUT), SPCOUT(nspchnl,N_SPCOUT)
      double precision reflfrac, w, ho_rt, ss

      integer IC2(N_IC), jwl, jcl, i

!     ss is a very small number subtracted so that leaf abedo would not be unity
!        which causes instability in the twostream model implemented in layer.f      
      ss = 1.0e-9
!
!     copy data
!      write(*,*)'copying data'
      do jcl = 1,nclmax
        do i = 1,N_XC
          XC2(jcl,i)=XC(jcl,i)
        enddo
      enddo
      do jwl = 1,nspchnl
        do i = 1,N_SPCIN
          SPCIN2(jwl,i)=SPCIN(jwl,i)
        enddo
!       initialize correction factor to a very wrong number
        cf(jwl) = -1.0e7
      enddo
      do i = 1,N_IC
        IC2(i)=IC(i)
      enddo
!      write(*,*)'corrfact: setting up forest data'
!     modify the needed values
      do jcl = 1, IC(1)
!       eliminate specular reflectance
        XC2(jcl,12) = 0
        do jwl = IC(10), IC(11)
!         leaf albedo
          w = (SPCIN(jwl,7+jcl)+SPCIN(jwl,7+nclmax+jcl))
          if (w .gt. 0 ) then
!           fraction of weighted reflectance to weighted total scattering, F = refl/(refl+tran)
!            set leaf R=F and T=1-F, so leaf W=1
            reflfrac = SPCIN(jwl,7+jcl) / w
!           make sure that leaves are not completely non-absorbing 
!                 as this causes a singularity in layer.f
            SPCIN2(jwl,7+jcl) = reflfrac*(1-ss)
            SPCIN2(jwl,7+nclmax+jcl) = (1.0-reflfrac)*(1-ss)
          else
!           black leaves, compute correction factor anyway
            SPCIN2(jwl,7+jcl) = 0.5*(1-ss)
            SPCIN2(jwl,7+nclmax+jcl) = 0.5*(1-ss)
          endif
!         set branch and trunk R to 1
          SPCIN2(jwl,7+3*nclmax+jcl) = 1.0*(1-ss)
          SPCIN2(jwl,7+4*nclmax+jcl) = 1.0*(1-ss)
!         in the future, also leaf abaxial reflectance needs to be set
!         SPCIN2(jwl,7+2*nclmax+jcl) = (1.0-reflfrac)*(1-ss)
        enddo
      enddo
!     set black ground and direct incidence; set cf to unity (no correction)
      do jwl = IC(10), IC(11)
        SPCIN2(jwl,2) = 0.0
        SPCIN2(jwl,3) = 0.0
        SPCIN2(jwl,4) = 0.0
        SPCIN2(jwl,5) = 0.0
        SPCIN2(jwl,6) = 1.0
        SPCIN2(jwl,8+5*nclmax) = 1.0
      enddo
!     set the flags to recalculate everything as we cannot assume anything has been computed
      IC2(4) = 1
      IC2(5) = 1
      IC2(7) = 1
      IC2(8) = 1
      IC2(9) = 1
c     set ijob to flux computation
      IC2(12+nclmax) = -1
      call comprt(XC2, IC2, SPCIN2, XOUT, SPCOUT)

!     calculate the correction factor 
      do jwl = IC(10), IC(11)
        ho_rt =  SPCOUT(jwl,5) + (SPCOUT(jwl,2)-SPCOUT(jwl,7))
        cf(jwl)=( 1-(SPCOUT(jwl,3)+SPCOUT(jwl,7)+XOUT(18)) ) / ho_rt
        write(*,'(a,f6.1,a,f10.6)') 'wl (nm)=',SPCIN(jwl,1),
     &            ' correction_factor=', cf(jwl)
      write(*,'(5f6.3)')SPCOUT(jwl,5),(SPCOUT(jwl,2)-SPCOUT(jwl,7)),
     &   SPCOUT(jwl,3),SPCOUT(jwl,7),XOUT(18) ! XXX
!    XXX print leaf R and T , and effective element R and T
      write(*,'(8f7.3)') SPCIN2(jwl,7+1), SPCIN2(jwl,7+nclmax+1),
     & SPCIN2(jwl,7+2), SPCIN2(jwl,7+nclmax+2),
     & SPCIN2(jwl,7+3), SPCIN2(jwl,7+nclmax+3),
     & SPCOUT(jwl,8),SPCOUT(jwl,9)
      enddo
      return
      end