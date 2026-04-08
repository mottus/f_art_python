!     part of the frt distribution, only used by python 
!     function to read the new input file format from python
!     depends on rd_cfm.f
      subroutine xd_cfm( fname, XC, IC, SPCIN, DESC, lerr )
!F2PY INTENT(IN) :: fname
!F2PY INTENT(OUT) :: XC, IC, SPCIN, DESC, lerr
compile: f2py.exe -c --compiler=mingw32 -m xd_cfm rd_cfm.o xd_cfm.f
c     A wrapper for rd_cfm for opening the input file
c     The file is usually opened in frt.f; this function
c     allows to call rd_cfm from python
c     IN: fname: file name
c     OUT: lerr: (logical) whether an error was detected
      implicit none
      include 'cf_new.h'
      integer N_XC, N_IC, N_SPCIN, iwl, jcl, j
      parameter ( N_XC = 30 )
      parameter ( N_IC = 20+nclmax )
      parameter ( N_SPCIN = 8+5*nclmax )
      integer IC(N_IC)
      double precision XC( nclmax, N_XC )
      double precision SPCIN( nspchnl, N_SPCIN )
c     DESC: stand name and class species' names
      character*77 DESC(nclmax+1)

      character*255 fname
      logical lerr
      integer ire
c     make the parameter in frtpar.h available in python via a /frtpar/
c      frtpar.h is read in cf_new.h
c     different names have to be used for parameters, common block variables and dummy variables
      integer fnclmax,fncub,fnspchnl,fnknotm,fnphim
      common /frtpar/ fnclmax, fncub, fnspchnl, fnknotm, fnphim
      
c     copy parameters from frtpar.h to /frtpar/      
      fnclmax = nclmax
      fncub = ncub
      fnspchnl = nspchnl
      fnknotm = nknotm
      fnphim = nphim

c     initialize outputs to zero
      do iwl=1,nspchnl
        do j = 1, N_SPCIN
          SPCIN(iwl,j) = 0
        enddo
      enddo
      do jcl = 1, nclmax
        DESC(jcl) = ' '
        do j=1,N_XC
          XC(jcl,j)=0.0
          XC(1,1)=0.0
        enddo
      enddo
      DESC(nclmax+1) = ' '
      do j=1,N_IC
        IC(j) = 0
      enddo

      ire = 12
      wl_set = .false.
      write(*,*)"xd_cfm opening file ", trim(fname)
      open (ire, file = fname, status = 'old')
      call rd_cfm(ire, XC, IC, SPCIN, DESC, lerr)
      close(ire)
      write(*,*)"xd_cfm finished reading configuration file."
      return
	  ire = 13
	  write(*,*) "jjj"
      end  ! subroutine xd_cfm( )
