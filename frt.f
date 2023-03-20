      program frt
c   Forest reflectance and transmittance
c   Based on the original fortran code of
c   (C) Andres Kuusk & Tiit Nilson
c   Kuusk, A. and Nilson, T. A directional multispectral forest
c   reflectance model, Remote Sens. Environ., 2000, 72(2):
c   244-252.
c
      implicit none
      include 'frtpar.h'
c
      integer N_XC, N_IC, N_SPCIN, N_XOUT, N_SPCOUT
      parameter ( N_XC = 30 )
      parameter ( N_IC = 20+nclmax )
      parameter ( N_SPCIN = 8+5*nclmax )
      parameter ( N_XOUT = 20 )
      parameter ( N_SPCOUT = 9 )
            
      character*77 clfile, oufile, fname
      character*260 buff ! temp. buffer for determining inp.file type
      character*77 DESC(nclmax+1) !stand name and class species' names
      real ctime, ct_arr(2) ! used for calling ETime at the end of program

      integer iwr, iwl, ivz, ire, ier, itemp, ijob
      integer j, jcl
      integer icorr ! wl index for printout of correction factor
      
      integer IC(N_IC), IC_temp(N_IC)
      double precision XC( nclmax, N_XC )
      double precision SPCIN( nspchnl, N_SPCIN )
      double precision XOUT( N_XOUT )
      double precision SPCOUT( nspchnl, N_SPCOUT )
      double precision hcrown, cf(nspchnl)
      double precision vz0 ! the view zenith angle range for ijob=2
      double precision azC ! the complementary view azimuth angle for ijob=2
      logical ltemp

!
!     XC, IC and SPCIN are fully described in comprt.f
!     IC has some additional elements not used in comprt, e.g.
! IC(2): correction wavelength index:
!           -1=calculate correction factor for each wavelength
!            0=no correction
!      integer=wavelength index
! IC(12+nclmax): ijob

c   =================================================================
c   Initialize variables
      pi = acos(-1.d0)
      dr = pi/180.d0
c     actually, pi and dr are also filled in comprt
      iwr    = 1 ! output file handle
      ire    = 2 ! input file handle
      ier    = 0
      l_volint = .false.
c     set outputs to zero. Inputs are zet to 0 in rd_cfm
      do iwl=1,nspchnl
        do j = 1,N_SPCOUT
          SPCOUT(iwl,j) = 0
        enddo
      enddo
      do j = 1, N_XOUT
        XOUT(j) = 0
      enddo

c   =================================================================
c   Read IJOB and forest parameters
c
      call getarg(1, clfile)
      call getarg(2, oufile)
      if (clfile .eq. ' ') then
         write (0, '(a)')  ' Name of input file '
         read (*, '(a)')  clfile
      endif
      if (oufile .eq. ' ') then
         write (0, '(a)')  ' Name of output file '
         read (*, '(a)') oufile
      endif
c
      ire = 2
      fname = clfile
      open (ire, file = fname, status = 'old', err = 99)
      read (ire, '(A)', err = 98) buff
      if ( buff(1:15) .ne. 'FRT CONFIG FILE' ) then
        write (*, '(a)')  'Not FRT config file in new format: ', fname
        stop 'FRT'
      endif

      call fndnb( fname, itemp )
      write(*,'(2a)') 'Reading input file (new format) ', fname
      rewind(ire)
      call rd_cfm( ire, XC, IC, SPCIN, DESC, ltemp )
      close (ire)
      if ( ltemp ) goto 98 ! error has occurred in rd_cfm
      ijob = IC(12+nclmax)
      write(*,*) 'Input read, job ID=',ijob

c     Earlier, subroutines used the following t_hknot: 
c     1) view zenith angle; 2 -- nknot: quadrature;
c     nknotm+1 -- nknotm+n_sun: sun zenith angles (filled below)
c     thknot is not used anymore. If some such old routines are 
c     bought back to life, the following code needs to be executed
c      nknot = nquad_t+1
C       nphi = nquad_p
C       t_hknot(1) = 0.d0
C       p_14(1) = 0.d0
C       p_phi = 2*pi / nphi
C       (t_hknot(i+1)=xquad_p(i),i,1,nquad_t)
C       (p_14(i+1)=wght_q(i)*nquad_p/2,i=1,nquad_t)
C       (p_hknot(i)=xquad_p(i),i=1,nquad_p)
      
c
c   =================================================================
c   CHECKS
c
c   Check crown length
      do jcl = 1, IC(1)
        if (IC(11+jcl). gt. 0) then
c         are crowns ellipsoids (l_elli)?
          if (XC(jcl, 3) .gt. XC(jcl, 2)) then
            XC(jcl, 3) = XC(jcl, 2)
            print *, 'jcl = ', jcl
            print *,'Warning: tree height is less than crown length,'
            print *, 'the crown length is adjusted.'
          endif
        else
          hcrown = XC(jcl, 3) + XC(jcl, 4)
          if (XC(jcl, 2) .lt. XC(jcl, 3)) then
            XC(jcl, 3) = XC(jcl, 2)
            XC(jcl, 4) = 0.d0
            print *, 'jcl = ', jcl
            print *, 'Warning: tree height is less than crown length,'
            print *, ' -> crown length is adjusted.'
          else if (hcrown .gt. XC(jcl, 2)) then
            XC(jcl, 4) = XC(jcl, 2) - XC(jcl, 3)
            print *, 'tree class = ', jcl
            print *, 'Warning: tree height is less than crown length,'
            print *, ' -> crown length is adjusted.'
          endif
        endif
      enddo

c  =================================================================
c     Calculate correction factors

      icorr = -1 ! indicating no correction
      if (IC(2) .eq. 0 ) then
c       calculate correction for each wavelength separately: slow!
        write(*,*) 'Correction for each wavelength separately: slow!'
        icorr = IC(10) ! the first wavelength
        call corrfact( XC, N_XC, IC, N_IC, SPCIN, N_SPCIN, 
     & N_XOUT, N_SPCOUT, SPCIN(1,8+5*nclmax) )
        write(*,*) 'Correction factors calculated.'
      elseif ( IC(2) .gt. IC(11) ) then
c       the index is out of range
        write(*,*) 'No correction used'
c       set correction factors to 1
        do iwl = IC(10),IC(11)
          SPCIN(iwl,8+5*nclmax) = 1.0
        enddo
      else
c        the wavelength index used for correction is stored in IC(2)
        write(*,*) 'Correction factors for wl #',IC(2)
        do j=1,N_IC
          IC_temp(j) = IC(j)
        enddo
        IC_temp(10) = IC(2)
        IC_temp(11) = IC(2)
        icorr = IC(2)
        call corrfact( XC, N_XC, IC_temp, N_IC, SPCIN, N_SPCIN, 
     & N_XOUT, N_SPCOUT, cf )
        write(*,*) 'Correction factors calculated.'
        do iwl = IC(10),IC(11)
          SPCIN(iwl,8+5*nclmax) = cf( IC(2) )
        enddo
      endif
c

c  =================================================================
c
c
c       CALL comprt -- THAT'S WHERE THE REAL CODE IS
c
c  set flags to recompute everything
      IC(4) = 1
      IC(5) = 1
      IC(6) = 1
      IC(7) = 1
      IC(8) = 1
      IC(9) = 1

c  create output file and store stand data
      write(*,'(2a)') 'Opening output file ', oufile
c      open (iwr, file = oufile, status = 'unknown')
      open (iwr, file = oufile, status = 'replace')
      rewind iwr
      write (iwr,'(a1,a11,a)') '#','name ', DESC(1)
      write (iwr,'(a1,a11,i4)') '#','age', nint(XC(1,25))
      
      if ( ijob .lt. 2 ) then   ! ------------------------------------
        call comprt(XC, IC, SPCIN, XOUT, SPCOUT)
!       single model run
!       store stand structural data in output file 
        write (iwr,'(a1,a11,F7.3)') '#','LAIeff_sh', XOUT(10)
        write (iwr,'(a1,a11,F7.3)') '#','LAI', XOUT(11)
        write (iwr,'(a1,a11,F7.3)') '#','PAI', XOUT(11)+XOUT(9)
        write (iwr,'(a1,a11,F7.3)') '#','crown_cl', XOUT(12)
        write (iwr,'(a1,a11,F7.3)') '#','can_cover', XOUT(13)
        write (iwr,'(a1,a11,F7.3)') '#','DIFN', XOUT(14)
        write (iwr,'(a1,a11,F7.3)') '#','PAI_opt', XOUT(15)
        write (iwr,'(a1,a11,F7.3)') '#','p_Pola', XOUT(16)
        write (iwr,'(a1,a11,F7.3)') '#','gap_Sun', XOUT(18)
        write (iwr,'(a1,a11,F7.3)') '#','LAI40deg', XOUT(20)
        write (iwr,'(a1,a11,F7.3)') '#','gap_view', XOUT(17)
        write (iwr,'(a1,a11,F7.3)') '#','gap_bidi', XOUT(19)
        if ( icorr .ne. -1 ) then
          write (iwr,'(a1,a7,I0,F7.3)') '#','Corr@', 
     &      int(SPCIN(icorr,1)), SPCIN(IC(10),8+5*nclmax)
        endif
!       write header line to output file
        write (iwr, '(a1, a5, 7(a12), 2(a8))')
     &    '#','wl', 'HDRF', 'HDTF', 
     &    'R1_cr' , 'R1_gnd' , 'Rhi_cr' , 'Rhi_gnd', 'T1',
     &    'R_elem', 'T_elem'
        do iwl = IC(10),IC(11)
          write (iwr,'(f6.1, 7(E12.4),2(F8.3))') SPCIN(iwl,1),
     &     ( SPCOUT(iwl,j), j=1,9)
      enddo
      
      elseif (ijob .eq. 2) then ! ------------------------------------
!       single model run is not enough, loop over view angles
        vz0 = XC(1,21)
        azC = XC(1,22) + pi ! the complement of view azimuth
        if (azC .gt. 2*pi) azC=azC-2*pi ! make sure to be < 2pi

        do ivz = 0, floor(2*vz0)
          XC(1,21) = vz0-ivz*vz0
          if ( XC(1,21) .lt. 0 ) then
            XC(1,21) = -XC(1,21)
            XC(1,22) = azC
          endif
          call comprt(XC, IC, SPCIN, XOUT, SPCOUT)
          if ( ivz .eq. 0 ) then
!           first view angle, store all stand data,
!           incl. view angle-independent data, in output file 
            write (iwr,'(a1,a11,F7.3)') '#','LAIeff_sh', XOUT(10)
            write (iwr,'(a1,a11,F7.3)') '#','LAI', XOUT(11)
            write (iwr,'(a1,a11,F7.3)') '#','PAI', XOUT(11)+XOUT(9)
            write (iwr,'(a1,a11,F7.3)') '#','crown_cl', XOUT(12)
            write (iwr,'(a1,a11,F7.3)') '#','can_cover', XOUT(13)
            write (iwr,'(a1,a11,F7.3)') '#','DIFN', XOUT(14)
            write (iwr,'(a1,a11,F7.3)') '#','PAI_opt', XOUT(15)
            write (iwr,'(a1,a11,F7.3)') '#','p_Pola', XOUT(16)
            write (iwr,'(a1,a11,F7.3)') '#','gap_Sun', XOUT(18)
            write (iwr,'(a1,a11,F7.3)') '#','LAI40deg', XOUT(20)
            write (iwr,'(a1,a11,F7.3)') '#','gap_view', XOUT(17)
            write (iwr,'(a1,a11,F7.3)') '#','gap_bidi', XOUT(19)
            if ( icorr .ne. -1 ) then
              write (iwr,'(a1,a7,I0,F7.3)') '#','Corr@', 
     &          int(SPCIN(icorr,1)), SPCIN(IC(10),8+5*nclmax)
            endif
!           write output file header line
            write (iwr, '(a1, a6, 2(a7), 7(a11), 2(a7))')
     &        '#', 'polar', 'azim', 'wl', 'HDRF', 'HDTF', 
     &        'R1_cr' , 'R1_gnd' , 'Rhi_cr' , 'Rhi_gnd', 'T1',
     &        'R_elem', 'T_elem'
          else
!           write only angle-dependent structural data
            write (iwr,'(a1,a11,F7.3)') '#','gap_view', XOUT(17)
            write (iwr,'(a1,a11,F7.3)') '#','gap_bidi', XOUT(19)
          endif
!         save simulation results
          do iwl = IC(10),IC(11)
            write (iwr,'(3(F7.1) 7(E11.4),2(F7.3))')
     &        XC(1,21)/dr, XC(1,22)/dr, SPCIN(iwl,1),
     &        ( SPCOUT(iwl,j), j=1,9 )   
          enddo
        enddo
      else 
        write(*,*) 'FRT: unknown ijob', ijob
      endif
c  =================================================================
!     close output file
      write (iwr,*) ''
      close(iwr) 

      call ETime( ct_arr, ctime )
      write(*,'(3(a,f8.1))') 'Time spent: ', ctime,
     &   '; of that user: ',ct_arr(1),' system:', ct_arr(2)
      print *, char(7) ! beep
      stop 'FRT'

c  End program with an error message

98    write (*, '(3a)')
     &      'Error in file ', fname
      stop 'FRT'
99    write (*, '(a)')  'File not found: ', fname
      stop 'FRT'

      end
