c     functions to read the new input file format
c     part of FRT
c
      subroutine xd_cfm( fname, XC, IC, SPCIN, DESC, lerr )
c     A wrapper for rd_cfm for opening the input file
c     The file is usually opened in frt.f; this function
c     allows to call rd_cfm from python
c     IN: fname: file name
c     OUT: lerr: (logical) whether an error was detected
CF2PY INTENT(IN) :: fname
CF2PY INTENT(OUT) :: XC, IC, SPCIN, DESC, lerr
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
      open (ire, file = fname, status = 'old', err = 99)
      call rd_cfm(ire, XC, IC, SPCIN, DESC, lerr)
      close(ire)
      write(*,*)"xd_cfm finished reading configuration file."
      return 
99    write (*, '(a)')  'File not found: ', trim(fname)
      return 
      end  ! subroutine xd_cfm( )

      subroutine rd_cfm( ifh, XC, IC, SPCIN, DESC, lerr )
c
CF2PY INTENT(IN) :: ifh
CF2PY INTENT(OUT) :: lerr, XC, IC, SPCIN, DESC
c     IN: ifh: input file handle
c     XC, IC, SPCIN model input parameter matrices
c     DESC: character array for stand description and class species' names
c     lerr: error flag
c     subroutine to read the configuration file

      implicit none
      include 'cf_new.h'      
      integer ifh
c     ifh: input file handle
      double precision XC(nclmax,*), SPCIN(nspchnl,*)
      integer IC(*)
      character*(*) DESC(*)
      logical lerr
c     lerr: (logical) whether an error was detected
      
      character*260 buff ! buffer for holding a line of input file
      character*32 sectyp, secnam, secpar
      integer i, nb, ir, itemp1, itemp2, itemp3, ierr
      logical lsec  
c       whether we're in the middle of reading a section
      logical lu 
c      whether sth useful was found in the loop
      logical lo 
c      whether there is still some input data missing
      logical loop
c
c     reqseq defined in cf_new.h.
c     Columns: type | name | the name of the SECTION requiring it
      do i = 1,nsecm
        reqsec(i,1)=' '
        reqsec(i,2)=' '
        reqsec(i,3)=' '
      enddo
      reqsec(1,1)='frt'
      lo = .true.
      lerr = .false.
      cf_err = .false.
      cf_warn = .false.
      wl_set = .false.
      wldelta = 0.25
c     when reading spectral properties from files, the maximum allowed
c            difference between model wl and file wl is wldelta [nm]
c     cln: current line number
      cln = 0 
      lsec = .false. 
c     the current tree class has to be initialized to zero, increased in rd_ell
      nn = 0
c
c     ----------------------------------------------------------
c     loop to reread input file until all required data has been
c     read or we've found that there's definitely sth missing.
c     the program does not buffer read SECTORs, so if a SECTOR
c     locateed closer to the beginning of the file is referred to,
c     the whole file has to be reread
c     Whether a SECTION is required or not is looked up from
c     reqsec; reqsec is filled as SECTIONs are being read
c
c     === OUTER UNTIL-LOOP TO READ THE FILE MANY TIMES OVER ===
 100  continue
      lu = .false.
c     ----------------------------------------------------------
c     inner loop to read the whole input file
      loop = .true.
      do while (loop)
        read (ifh, '(A)', IOSTAT=ierr) buff
        cln = cln + 1
        if (ierr .eq. 0) then
c         delete comments from the end of line & leading spaces
c           and report first non-blank character
          call fndnb( buff, nb )
c         we're looking for SECTIONs...
          if ( nb .ge. 1 ) then 
c           this is not a blank line
            if ( buff(1:7) .eq. 'SECTION ' ) then
              lsec = .true.
c             check whether we're interested in this SECTION
              call spliln( buff, sectyp, secnam) ! split
              call isreq( sectyp, secnam, ir ) ! do we need it?
              if ( ir .gt. 0 )  then 
c               this section has to be read
                secpar = reqsec(ir,3)
c               secpar contains the name of the parent section requesting this one
                call rdsec( ifh, sectyp, secpar, XC, IC, SPCIN, DESC ) ! read
                call remreq( ir ) ! remove request for the read section
                lu = .true.
                lsec = .false.
              endif
            elseif ( buff(1:4) .eq. 'END' .and. .not. lsec ) then
                write(*,*) cln,
     &            ' WARNING: END although not in a SECTION'
                lsec = .false.
            endif
c           I/O errors and warnings
            if ( cf_warn .or. cf_err ) then
              call fndlnb( sectyp, itemp1 )
              call fndlnb( secnam, itemp2 )
              if ( cf_warn ) then
                write(*,*) '  warning in SECTION ',
     &            sectyp(1:itemp1),' ',secnam(1:itemp2)
                cf_warn = .false.
              endif
              if ( cf_err ) then
                write(*,*) '  error in SECTION ',
     &            sectyp(1:itemp1),' ',secnam(1:itemp2)
                cf_err = .false.
                lerr = .true.
              endif
            endif
          endif 
        else 
c         that is, (ierr .ne. 0), exit the loop due to EOF
          loop = .false.
        endif
      enddo
c     ----------------------------------------------------------
c     end of while loop, we are here only when reading produces an error at EOF
c
c     check whether there's sth that was not found
      lo = .false.
      do i = 1,nsecm
        if ( reqsec(i,1) .ne. ' ' ) then
          lo = .true.
        endif
      enddo
c     do we have to reread the file?
      if ( lo ) then
        if ( lu ) then
          rewind ifh          
          cln = 0 ! reset line number
          goto 100
        else
          write(*,*) 'ERROR, data missing:'
          do i = 1,nsecm
            if ( reqsec(i,1) .ne. ' ' ) then
              call fndlnb( reqsec(i,1), itemp1 )
              call fndlnb( reqsec(i,2), itemp2 )
              call fndlnb( reqsec(i,3), itemp3 )
              write(*,*) ' SECTION ', reqsec(i,1)(1:itemp1),
     &          ' ', reqsec(i,2)(1:itemp2),
     &          '   <- SECTION ', reqsec(i,3)(1:itemp3)
              lerr = .true.
            endif
          enddo
        endif
      endif
c     === END OF OUTER UNTIL-LOOP ================================
c        no need to re-read the file, we have all we need (or an error)
      write(*,'(a,f6.0,a,f6.0,a)') 'Usable wavelength range: ', 
     &  SPCIN(IC(10),1),' -- ', SPCIN(IC(11),1), ' nm'
c     check if correction wavelength index, IC(2) is inside the range
c     assume that it refers to a wavelength counted from 
c        the beginning of the wavelength array, SPCIN(1,1)
      if ( IC(2) .gt. 0 ) then
        if ( IC(2) .lt. IC(10) .or. IC(2) .gt. IC(11) ) then
          write(*,*) 'WARNING: correction_wl out of range,'//
     &      ' NO CORRECTION'
c         set the index to beyond the existing range
          IC(2) = IC(11)+100
        endif
      endif
c     a final touch: the p_khair (or pkhair) parameter is not read from a file, and needs to be set to 1
      XC(1,24) = 0.1 ! unclear, why this needs to be the value, it affects specular gamma-function in gmfres
      return
      end  ! subroutine rd_cfm
c
c
      subroutine spliln( buff, sectyp, secnam )
c  splits the input buffer into words (finds SECTION type and name)
c  IN: buff: input buffer declared in cf_new.h
c  OUT: sectyp, secname: type and name of the section
      implicit none
      include 'cf_new.h'
      character*(*) buff,sectyp, secnam
      
      character*260 buff2
      integer i
c
c     setting outputs to blank: if nothing is found in buff, we don't want the old value anyway
      sectyp=' '
      secnam=' '
c
c     assuming that buff(1:8)='SECTION '
      buff2 = buff(9:len(buff))
c     1st word
      call fndnb( buff2, i ) ! delete leading blanks & comments
      if ( i .le. 0 ) goto 99
      i = index( buff2, ' ' )
      if ( i .le. 0  ) goto 99
c       There should always be blanks at the end of buff2 due to
c       chopping off its beginning.
      sectyp = buff2(1:i-1)
      buff2 = buff2(i+1:len(buff2))
c     2nd word
      call fndnb( buff2, i ) ! delete leading blanks
      if ( i .le. 0 ) goto 99
      i = index( buff2, ' ' )
      if ( i .eq. 0 ) i = len(buff2) ! this should never happen
      secnam = buff2(1:i-1)
 99   continue
      return
      end  ! subroutine spliln
c
c
      subroutine splil2( buff, sectyp, secnam )
c  gets the second and third words from buff
c  IN: buff: input buffer
c  OUT: sectyp, secname: type and name of the section,
c    2nd & 3rd word of line, respectively. splil2 could be used
c    also instead of spliln, but spliln was written earlier
c
      implicit none
      character*(*) buff, sectyp, secnam
      
      character*260 buff2
      integer i
c
c     setting outputs to blank: if nothing is found in buff, we don't want the old value anyway
      sectyp=' '
      secnam=' '
c
c     assuming that blanks have been stripped
      i = index( buff, ' ' )
      if ( i .le. 0  ) goto 99 ! no blanks
      buff2 = buff(i:len(buff)) ! 2nd word
      call fndnb( buff2, i ) ! delete leading blanks
      if ( i .le. 0 ) goto 99
      i = index( buff2, ' ' )
      if ( i .le. 0  ) goto 99
c       There should always be blanks at the end of buff2 due to
c       chopping off its beginning.
      sectyp = buff2(1:i-1)
      buff2 = buff2(i+1:len(buff2))!  2nd word
      call fndnb( buff2, i ) ! delete leading blanks
      if ( i .le. 0 ) goto 99
      i = index( buff2, ' ' )
      if ( i .eq. 0 ) i = len(buff2) ! this should never happen
      secnam = buff2(1:i-1)
 99   continue
      return
      end  ! subroutine splil2
c
c
      subroutine fndnb( buff, nb )
c
c     removes leading spaces & trailing comments and gives the
c     location of the first non-blank character of the original
c     unmodified string buff
c  IN/OUT:  buff: input string (comments are removed)
c  OUT: nb: location of 1st nonblank (-1 if not found)
      implicit none
      character*(*) buff
      integer nb, nc
c
      nb = 1
 110  continue 
c     loop to find first non-blank position
        if ( buff(nb:nb) .eq. ' ' .OR. buff(nb:nb) .eq. char(9) ) then
c          char(9) = TAB
          if ( nb .lt. len(buff) ) then
            nb = nb + 1
            goto 110
          else
            nb = -1
          endif
        elseif ( buff(nb:nb) .eq. '#' ) then     
c         treat as a blank line
          buff = ' ' 
          nb = -1
        endif
c     ----- end loop
      if ( nb .gt. 0 ) then 
c       if not empty, check for comments
        nc = nb+1 ! nc: location of first '#'
 120    continue 
c       loop to check for comments
          if ( buff(nc:nc) .ne. '#' ) then
            nc = nc + 1
            if ( nc .le. len(buff) ) goto 120
          endif
          buff = buff(nb:nc-1)
c       ----- end loop
      endif
c
      return
      end  ! subroutine fndnb
c
c
      subroutine rdsec( ifh, sectyp, secpar, XC, IC, SPCIN, DESC )
c   call the correct subroutine given section at hand
c   IN: ifh: input file handle
c    sectyp the type of section at hand
c    secpar: the name of the SECTION requiring it
c        secpar is not currently used (2022)
c    XC, IC, SPCIN model input parameter matrices
c    DESC: character array for stand description and class species' names
c      switch between subroutines depending on SECTION type
      implicit none
      include 'cf_new.h'
      integer ifh
      character*(*) sectyp, secpar
      double precision XC(nclmax,*), SPCIN(nspchnl,*)
      integer IC(*)
      character*(*) DESC(*)
c
      if ( sectyp .eq. 'frt' ) then
        if ( wl_set ) then
c          read the rest of the parameters, except the used wavelengths        
          call rd_frt2( ifh, XC, IC, SPCIN, DESC )
        else
c         first time around, read the used wavelengths
          call rd_frt1( ifh, XC, IC, SPCIN )
        endif
      elseif ( sectyp .eq. 'ellipsoid' ) then
        call rd_ell( ifh, 1, XC, IC, SPCIN, DESC )
      elseif ( sectyp .eq. 'cone' ) then
        call rd_ell( ifh, 2, XC, IC, SPCIN, DESC )
      elseif ( sectyp .eq. 'difsky' ) then
        call rd_sky( ifh, IC, SPCIN  )
      elseif ( sectyp .eq. 'geometry' ) then
        call rd_geo( ifh, XC, IC )
      elseif ( sectyp .eq. 'spectrum' ) then
        call rd_spc( ifh, IC, SPCIN )
      elseif ( sectyp .eq. 'discrete') then
        call rd_dsc( ifh, IC, SPCIN )
      elseif ( sectyp .eq. 'lambert' ) then
        call rd_lmb( ifh, IC, SPCIN )
      elseif ( sectyp .eq. 'GLcubature' ) then
        call rd_glc( ifh, IC )
      else
        write(*,*)' ERROR: dont know what to do with SECTION ',sectyp
        cf_err = .true.
      endif
      return
      end  ! subroutine rdsec
c
c
      subroutine rd_glc( ifh , IC )
c     the actual input if SECTION is GLcubature
c  IN: ifh: input file handle
c  OUT: IC: model input matrix
      implicit none
      include 'cf_new.h'
      integer IC(*), ifh
c
      character*260 buff
      integer nflds, ierr, i
      parameter (nflds=2)
      logical loop, ll, rv 
c   flags
      character*16 fields( nflds ) 
c   field names
      logical filled( nflds ) 
c   whether a value has been found in the file
      character*30 field
c
      data fields / 'theta', 'phi' /
      data filled / nflds * .false. /
c
c     -----------------------------------------------------------------
c     LOOP for reading lines from input file
      loop = .true.
      ierr = 0
      do while (loop)
        read (ifh, '(A)', IOSTAT=ierr) buff
        cln = cln + 1
c       error check
        if ( ierr .ne. 0 ) then
          write (*,*) cln,' WARNING: unexpected end of file'
c         exit the loop immediately
c         variables will be checked after the loop, so if all
c         important fields were filled, there is no reason to
c         create an error
          loop = .false.
          exit
        else 
c         strip blanks and check line
          call chsl( ifh, buff, ll, loop )
        endif
        if ( ll ) then
          read(buff,*) field
c         read the value into appropriate global variable
          do i=1,nflds
c           first data fields
            if ( fields(i) .eq. field ) then
              rv = .false.
              if ( i .eq. 1 ) then  ! nquad_t (theta)
                read(buff,*,iostat=ierr) field, IC(13+nclmax) 
                rv = .true.
              elseif ( i .eq. 2 ) then ! nquad_p (phi)
                read(buff,*,iostat=ierr) field, IC(14+nclmax) 
                rv = .true.
              endif
              if ( rv ) then  
c               error checking
                if ( ierr .ne. 0 ) then
                  write(*,*) cln,' Warning: Subroutine rd_glc '//
     &              ' could not read ', field
                  cf_warn = .true.
                else
                  filled(i) = .true.
                endif
              endif
c             break out of the "do i=1,nflds" loop
              exit
            endif
          enddo
          if (i .gt. nflds) then
            write(*,*) cln,' WARNING: unknown field: ', field
            cf_warn = .true.
          endif
        endif
      enddo
c     END  ----  do while (loop)
c     check whether there's something missing
      do i = 1, nflds
        if ( .not. filled(i) ) then
          write(*,*) 'ERROR: required field missing: ',fields(i)
          cf_err = .true.
        endif
      enddo
c
c     check
      if ( .not. cf_err ) then
        if ( IC(13+nclmax) .le. 0 .or. IC(14+nclmax) .lt. 0 ) then
          write(*,*) 'ERROR: less than 0 cubature knots'
          cf_err = .true.
        endif
        if ( (IC(13+nclmax) .gt. nknotm) .or. 
     &      (IC(14+nclmax) .gt. nphim) )  then
          write(*,*) "ERROR: too many cubature knots, 
     &         limited by nknotm and nphim in frtpar.h"
          cf_err = .true.
        endif
      endif
c
      return
      end ! subroutine rd_glc
c
c
      subroutine rd_dsc( ifh, IC, SPCIN )
c     the actual input if SECTION is discrete
c  IN: ifh: input file handle; 
c  INOUT    IC, SPCIN: model input parameter vectors
      implicit none
      include 'cf_new.h'
      double precision SPCIN(nspchnl,*)
      integer IC(*), ifh
c
      character*260 buff
      logical loop, ll 
c   flags
      integer nchnl, ierr
c     nchnl: the number of channels already read
      nchnl = 0

c     LOOP for reading lines from input file -------------
      loop = .true.
      ierr = 0
      do while (loop)
        read (ifh, '(A)', IOSTAT=ierr) buff
        cln = cln + 1
c       error check
        if ( ierr .ne. 0 ) then
          write (*,*) cln,' WARNING: unexpected end of file'
c         exit the loop immediately
c         variables will be checked after the loop, so if all
c         important fields were filled, there is no reason to
c         create an error
          loop = .false.
          exit
        else 
c         strip blanks and check line
          call chsl( ifh, buff, ll, loop )
        endif
        if ( ll ) then
c         read the value into wavelength array, 1st column of SPCIN
          read(buff,*,iostat=ierr) SPCIN(nchnl+1,1)
          if ( ierr .ne. 0 ) then
            write(*,*) cln,' Warning: Subroutine rd_dsc'//
     &              'could not read wavelengths'
            cf_warn = .true.
          elseif (nchnl .eq. nspchnl ) then
c           maximum array length reached
            write(*,*) cln,' Warning: Subroutine rd_dsc'//
     &              'maximum number of wavelengths reached, stopping.'
            loop = .false.
          else
            nchnl = nchnl + 1
          endif
        endif
      enddo
c     END  ----  do while (loop)
c     check whether there's something missing
      IC(10) = 1
      IC(11) = nchnl
      if ( nchnl .eq. 0 ) then
        write(*,*) 'ERROR: no channels found'
        cf_err = .true.
      endif
      return
      end  ! subroutine rd_dsc
c
c
      subroutine rd_lmb( ifh, IC, SPCIN )
c
c     the actual input if SECTION is lambert (for ground reflectance)
c  IN: ifh: input file handle; 
c  INOUT:    IC, SPCIN: matrix of input parameters and spectra, respectively
      implicit none
      include 'cf_new.h'
      integer ifh, IC(*)
      double precision SPCIN(nspchnl,*)
      integer nflds
c
      character*260 buff
      parameter (nflds=2)
c       flags
      logical loop, ll, rv 
c       field names
      character*16 fields( nflds ) 
c       whether a value has been found in the file
      logical filled( nflds ) 
c       whether a value is compulsory
      logical reqd( nflds ) 
      character*30 field
c      requested file name
      character*80 fname 
      double precision tmp_dbl
      integer i, i3, ierr
c
      data fields / 'reflectance' , 'file' /
      data filled / nflds * .false. /
      data reqd / nflds * .true. /
c
c     -----------------------------------------------------------------
c     LOOP for reading lines from input file
      loop = .true.
      ierr = 0
      do while (loop)
        read (ifh, '(A)', IOSTAT=ierr) buff
        cln = cln + 1
c       error check
        if ( ierr .ne. 0 ) then
          write (*,*) cln,' WARNING: unexpected end of file'
c         exit the loop immediately
c         variables will be checked after the loop, so if all
c         important fields were filled, there is no reason to
c         create an error
          loop = .false.
          exit
        else 
c         strip blanks and check line
          call chsl( ifh, buff, ll, loop )
        endif
        if ( ll ) then
          read(buff,*) field
c         read the value into appropriate global variable
          do i=1,nflds
c           first data fields
            if ( fields(i) .eq. field ) then
              rv = .false.
              if ( i .eq. 1 ) then ! reflectance 
                read(buff,*) field, tmp_dbl 
                do i3 = IC(10), IC(11)
!                 set all ground reflectance factors to the same value. This may change in the future  
                  SPCIN(i3,2) = tmp_dbl
                  SPCIN(i3,3) = tmp_dbl
                  SPCIN(i3,4) = tmp_dbl
                  SPCIN(i3,5) = tmp_dbl
                enddo
                rv = .true.
c               remove other requirements
                reqd(2) = .false.
              elseif ( i .eq. 2 ) then ! file
                read(buff,*) field, fname 
                call rfile(fname, SPCIN(1,1), IC(10), IC(11),
     &            SPCIN(1,2), 1, ierr)
                do i3 = IC(10), IC(11)
!                 set all ground reflectance factors to the same value. This may change in the future 
!                    e.g., by allowing the file to have more than one column 
                  SPCIN(i3,3) = SPCIN(i3,2)
                  SPCIN(i3,4) = SPCIN(i3,2)
                  SPCIN(i3,5) = SPCIN(i3,2)
                enddo
                rv = .true. 
c               remove other requirements
                reqd(1) = .false.
              endif
              if ( rv ) then  
c               error checking
                if ( ierr .ne. 0 ) then
                  write(*,*) cln,' Warning: Subroutine rd_lmb'//
     &              ' could not read ', field
                  cf_warn = .true.
                else
                  filled(i) = .true.
                endif
              endif
c             break out of the "do i=1,nflds" loop
              exit
            endif
          enddo
          if (i .gt. nflds) then
            write(*,*) cln,' WARNING: unknown field: ', field
            cf_warn = .true.
          endif
        endif
      enddo
c     END  ----  do while (loop)
c     check whether there's something missing
        do i = 1, nflds
          if ( reqd(i) .and. .not. filled(i) ) then
            write(*,*) 'ERROR: required field missing: ',fields(i)
            cf_err = .true.
          endif      
        enddo
      return
      end ! subroutine rd_lmb
c
c
      subroutine rd_spc( ifh, IC, SPCIN )
c
c     the actual input if SECTION is spectrum
c  IN: ifh: input file handle
c  INOUT:   IC, SPCIN: model input matrices
      implicit none
      include 'cf_new.h'
      integer ifh
      double precision SPCIN(nspchnl,*)
      integer IC(*)
c
      character*260 buff
      integer i, nflds, jwl, ierr
      parameter (nflds=3)
      logical loop, ll, rv 
c       flags
      character*16 fields( nflds ) 
c       field names
      logical filled( nflds ) 
c       whether a value has been found in the file
      character*30 field
      double precision dwl
c
      data fields / 'start', 'step', 'number' /
      data filled / nflds * .false. /
c
c     -----------------------------------------------------------------
c     LOOP for reading lines from input file
      loop = .true.
      ierr = 0
      do while (loop)
        read (ifh, '(A)', IOSTAT=ierr) buff
        cln = cln + 1
c       error check
        if ( ierr .ne. 0 ) then
          write (*,*) cln,' WARNING: unexpected end of file'
c         exit the loop immediately
c         variables will be checked after the loop, so if all
c         important fields were filled, there is no reason to
c         create an error
          loop = .false.
          exit
        else 
c         strip blanks and check line
          call chsl( ifh, buff, ll, loop )
        endif
        if ( ll ) then
          read(buff,*) field
c         read the value into appropriate global variable
          do i=1,nflds
c           first data fields
            if ( fields(i) .eq. field ) then
              rv = .false.
              if ( i .eq. 1 ) then ! start
                read(buff,*,iostat=ierr) field, SPCIN(1,1) 
                rv = .true.
              elseif ( i .eq. 2 ) then ! step
                read(buff,*,iostat=ierr) field, dwl 
                rv = .true.
              elseif ( i .eq. 3 ) then ! number
                read(buff,*,iostat=ierr) field, IC(11) 
                rv = .true.
              endif
              if ( rv ) then
c               error checking
                if ( ierr .ne. 0 ) then
                  write(*,*) cln,' Warning: Subroutine rd_spc'//
     &              ' could not read ', field
                  cf_warn = .true.
                else
                  filled(i) = .true.
                endif
              endif
c             break out of the "do i=1,nflds" loop
              exit
            endif
          enddo
          if (i .gt. nflds) then
            write(*,*) cln,' WARNING: unknown field: ', field
            cf_warn = .true.
          endif
        endif
      enddo
c     END  ----  do while (loop)
c     check whether there's something missing
      do i = 1, nflds
        if (.not. filled(i) ) then
          write(*,*) 'ERROR: required field missing: ',fields(i)
          cf_err = .true.
        endif
      enddo
c     fill SPCIN first column from this data (i.e., generate the spectrum)
c     fill also IC(10) and IC(11), start and end indices. Initially, start at 1.
      IC(10) = 1 
      if ( .not. cf_err ) then
c       SPCIN(1,1) contains the first wavelength 
        do jwl = 2, IC(11)
          SPCIN(jwl,1) = SPCIN(1,1) + (jwl-1) * dwl
        enddo
        wldelta = 0.3*dwl
        wl_set = .true.
        write(*,"(a,F6.1,a,F6.1)")"rd_spc read spectrum, wl range",
     &   SPCIN(1,1),"--",SPCIN(IC(11),1)
      endif
c
      return
      end ! subroutine rd_spc
c
c
      subroutine rd_geo( ifh, XC, IC )
c
c     the actual input if SECTION is geometry
c  IN: ifh: input file handle
c  OUT:  XC, model input parameter matrix
      implicit none
      include 'cf_new.h'
      double precision XC(nclmax,*)
      integer ifh, IC(*)
      
      character*260 buff
      integer i, nflds, ijob, ierr
      parameter (nflds=5)
c       flags
      logical loop, ll, rv 
c       field names
      character*16 fields( nflds ) 
c       whether a value has been found in the file
      logical filled( nflds ) 
c       whether a value is compulsory
      logical reqd( nflds ) 
      character*30 field
c
      data fields / 'sunzenith', 'viewnadir',
     &  'viewazimuth', 'viewincrement', 'sunzeniths' /
c   note: sunzeniths added for backward compatibility, same information as sunzenith
c     only one solar angle is allowed since 2021
      data filled / nflds * .false. /
      data reqd / nflds * .true. /
c
c     find the required parameters
      pi = acos(-1.d0)
      ijob = IC(12+nclmax)
      if ( ijob .eq. -1 ) then
        reqd(2) = .false.
        reqd(3) = .false.
        reqd(4) = .false.
      elseif ( ijob .eq. 0 ) then
        reqd(4) = .false.
      elseif ( ijob .eq. 1 ) then
        reqd(4) = .false.
      elseif ( ijob .eq. 2 ) then
        reqd(2) = .false.
      endif
      reqd(5) = .false.
c  sunzeniths not required, if found, same data as in sunzenith
c
c     -----------------------------------------------------------------
c     LOOP for reading lines from input file
      loop = .true.
      ierr = 0
      do while (loop)
        read (ifh, '(A)', IOSTAT=ierr) buff
        cln = cln + 1
c       error check
        if ( ierr .ne. 0 ) then
          write (*,*) cln,' WARNING: unexpected end of file'
c         exit the loop immediately
c         variables will be checked after the loop, so if all
c         important fields were filled, there is no reason to
c         create an error
          loop = .false.
          exit
        else 
c         strip blanks and check line
          call chsl( ifh, buff, ll, loop )
        endif
        if ( ll ) then
          read(buff,*) field
c         read the value into appropriate global variable
          do i=1,nflds
c           first data fields
            if ( fields(i) .eq. field ) then
              rv = .false.
              if ( i .eq. 1 ) then ! sunzenith
                read(buff,*,iostat=ierr) field, XC(1,23)
                XC(1,23) = XC(1,23)*pi/180.0
                rv = .true.
              elseif ( i .eq. 2 ) then ! viewnadir
                read(buff,*,iostat=ierr) field, XC(1,21)
                XC(1,21) = XC(1,21)*pi/180.0
                rv = .true.
              elseif ( i .eq. 3 ) then  ! viewazimuth
                read(buff,*,iostat=ierr) field, XC(1,22)
                XC(1,22)  = XC(1,22) * pi /180.0
                rv = .true.
              elseif ( i .eq. 4 ) then  ! viewincrement
                read(buff,*,iostat=ierr) field, XC(1,26)
                XC(1,26) = XC(1,26)*pi/180.0
                rv = .true.
              elseif ( i .eq. 5 ) then ! sunzeniths, for backward compatibility
                read(buff,*,iostat=ierr) field, XC(1,23)
                XC(1,23) = XC(1,23)*pi/180.0
                rv = .true.
                reqd(1) = .false.
              endif
              if ( rv ) then
c               error checking
                if ( ierr .ne. 0 ) then
                  write(*,*) cln,' Warning: Subroutine rd_geo'//
     &              ' could not read ', field
                  cf_warn = .true.
                else
                  filled(i) = .true.
                endif
              endif
c             break out of the "do i=1,nflds" loop
              exit
            endif
          enddo
          if (i .gt. nflds) then
            write(*,*) cln,' WARNING: unknown field: ', field
            cf_warn = .true.
          endif
        endif
      enddo
c     END  ----  do while (loop)
c     check whether there's something missing
      do i = 1, nflds
        if ( reqd(i) .and. .not. filled(i) ) then
          write(*,*) 'ERROR: required field missing: ',fields(i)
          cf_err = .true.
        endif
      enddo
c
      if ( XC(1,26) .lt. 0.01*pi/180.0 ) then
        write(*,*)
     &      "Warning: view increment too small, setting to 1 degree"
        XC(1,26) = 1.0*pi/180.0
        cf_warn = .true.
      endif

      return
      end ! subroutine rd_geo
c
c
      subroutine rd_sky( ifh, IC, SPCIN )
c
c     the actual input if SECTION is difsky
c  IN: ifh: input file handle
c  INOUT: IC, SPCIN: matrix of input parameter and spectra matrix, respectively
      implicit none
      include 'cf_new.h'
      integer ifh, IC(*)
      double precision SPCIN(nspchnl,*)

      character*260 buff
c       flags
      integer nflds
      parameter (nflds=2)
      logical loop, ll, rv 
c       field names
      character*16 fields( nflds ) 
c       whether a value has been found in the file
      character*80 fname 
      logical reqd( nflds ) 
      logical filled( nflds ) 
      character*30 field
      double precision tmp_dbl
      integer i, ierr, iw
c
      data fields / 'SQ_ratio', 'file' /
      data filled / nflds * .false. /
      data reqd / nflds * .true. / 
c
c     -----------------------------------------------------------------
c     LOOP for reading lines from input file
      loop = .true.
      ierr = 0
      do while (loop)
        read (ifh, '(A)', IOSTAT=ierr) buff
        cln = cln + 1
c       error check
        if ( ierr .ne. 0 ) then
          write (*,*) cln,' WARNING: unexpected end of file'
c         exit the loop immediately
c         variables will be checked after the loop, so if all
c         important fields were filled, there is no reason to
c         create an error
          loop = .false.
          exit
        else 
c         strip blanks and check line
          call chsl( ifh, buff, ll, loop )
        endif
        if ( ll ) then
          read(buff,*) field
c         read the value into appropriate global variable
          do i=1,nflds
c           first data fields
            if ( fields(i) .eq. field ) then
              rv = .false.
              if ( i .eq. 1 ) then
                read(buff,*,iostat=ierr) field, tmp_dbl ! SQ_ratio
                do iw=IC(10),IC(11)
                  SPCIN(iw,6) = tmp_dbl
                enddo
                rv = .true.
                reqd(2) = .false.
              elseif ( i .eq. 2 ) then ! file
                read(buff,*) field, fname 
c               
                call rfile(fname, SPCIN(1,1), IC(10), IC(11),
     &            SPCIN(1,6), 1, ierr)
                if ( ierr .eq. 0 ) then 
                  rv = .true. 
c                 remove other requirements
                  reqd(1) = .false.
                endif
              endif
              if ( rv ) then
c               error checking
                if ( ierr .ne. 0 ) then
                  write(*,*) cln,' Warning: Subroutine rd_sky'//
     &              ' could not read ', field
                  cf_warn = .true.
                else
                  filled(i) = .true.
                endif
              endif
c             break out of the "do i=1,nflds" loop
              exit
            endif
          enddo
          if (i .gt. nflds) then
            write(*,*) cln,' WARNING: unknown field: ', field
            cf_warn = .true.
          endif
        endif
      enddo
c     END  ----  do while (loop)
c     check whether there's something missing
      do i = 1, nflds
        if ( reqd(i) .and. .not. filled(i) ) then
          write(*,*) 'ERROR: required field missing: ',fields(i)
          cf_err = .true.
        endif
      enddo
c
      return
      end ! subroutine rd_sky
c
c
      subroutine rd_ell( ifh, itype, XC, IC, SPCIN, DESC )
c     the actual input if SECTION is ellipsoid or cylinder+cone
c  IN: ifh: input file handle
c    itype: 1 = ellipsoid, 2 = cylinder+cone
c  INOUT: XC, IC, SPCIN model input parameter matrices
c  OUT: XC: model input parameter matrices
c   DESC: character array for stand description and class species' names
      implicit none
      include 'cf_new.h'
      integer ifh, itype
      double precision XC(nclmax,*), SPCIN(nspchnl,*)
      integer IC(*)
      character*(*) DESC(*)

      character*260 buff
      integer nflds, i, i1, i5, ierr, iw
      parameter (nflds=18)
      logical loop, ll, rv 
c       flags
      character*16 fields( nflds ) 
c       field names
      logical filled( nflds ) 
c       whether a value has been found in the file
      logical use_p
c       use_p: whether the reflectance spectra are given for a needle and we need to scale
c       it to a highel level (e.g., shoot) using p calculated from ssc
      character*30 field, stemp, stemp2
      character*30 st, sn 
c       requested SECTION type & name
c       the number of the crown class being currently read
      double precision p_xarr(nspchnl,2) 
c       used to read leaf refl & trans from file
c       used to read trunk & branch refl from file
      double precision rtemp, ttemp, p, w

      data fields / 'density', 'height', 'crownlength', 'cyllength',
     &   'crownradius', 'trunkdiameter', 'dlw', 'slw',
     &   'bailai', 'tdp', 'ssc', 'shl',
     &   'br_refl', 'tr_refl', 'species', 'leafmodel',
     &   'wax_cf', 'scale_needle' /
c
c     set tree class number (defined in cf_new.h, initialized in rd_cfm to 0)
      nn = nn +1 
      IC(11+nn) = 2 - itype 
c     cyllength is a valid parameter only if cyl+cone model is used
      if (itype .eq. 2) then
        fields(4) = 'cyllength'
      else
        fields(4) = ' '
        XC(nn,4) = 0
      endif
      do i = 1,nflds
        filled(i) = .false.
      enddo
c      default value exists for use_p
      use_p = .false.
      filled(18) = .true.
c
c     -----------------------------------------------------------------
c     LOOP for reading lines from input file
      loop = .true.
      ierr = 0
      do while (loop)
        read (ifh, '(A)', IOSTAT=ierr) buff
        cln = cln + 1
c       error check
        if ( ierr .ne. 0 ) then
          write (*,*) cln,' WARNING: unexpected end of file'
c         exit the loop immediately
c         variables will be checked after the loop, so if all
c         important fields were filled, there is no reason to
c         create an error
          loop = .false.
          exit
        else 
c         strip blanks and check line
          call chsl( ifh, buff, ll, loop )
        endif
        if ( ll ) then
          read(buff,*) field
c         read the value into appropriate global variable
          do i=1,nflds
            rv = .false.
            if ( fields(i) .eq. field ) then
c
c             first requests to read a file
              call splil2( buff, st, sn )
              if ( st .eq. 'file' ) then
                if ( i .eq. 13 ) then ! br_refl
                  call rfile(sn, SPCIN(1,1), IC(10), IC(11),
     &              SPCIN(1,7+3*nclmax+nn), 1, ierr)
                  rv = .true.
                elseif ( i .eq. 14 ) then ! tr_refl
                  call rfile(sn, SPCIN(1,1), IC(10), IC(11),
     &              SPCIN(1,7+4*nclmax+nn), 1, ierr)
                  rv = .true.
                elseif ( i .eq. 16 ) then ! leafmodel
                  call rfile(sn, SPCIN(1,1), IC(10), IC(11),
     &              p_xarr, 2, ierr)
                  do i5 = 1, nspchnl
                    SPCIN(i5,7+nn) = p_xarr(i5,1)
                    SPCIN(i5,7+nclmax+nn) = p_xarr(i5,2)
                  enddo
                  rv = .true.
                endif
              else
c                now data fields
                if ( i .eq. 1 ) then ! density
                  read(buff,*,iostat=ierr) field, XC(nn,1)
                  rv = .true.
                elseif ( i .eq. 2 ) then ! height
                  read(buff,*,iostat=ierr) field, XC(nn,2)
                  rv = .true.
                elseif ( i .eq. 3 ) then ! crownlength
                  read(buff,*,iostat=ierr) field, XC(nn,3)
                  rv = .true.
                elseif ( i .eq. 4 ) then ! cyllength
                  read(buff,*,iostat=ierr) field, XC(nn,4)
                  rv = .true.
                elseif ( i .eq. 5 ) then ! crownradius
                  read(buff,*,iostat=ierr) field, XC(nn,5)
                  rv = .true.
                elseif ( i .eq. 6 ) then ! trunkdiameter
                  read(buff,*,iostat=ierr) field, XC(nn,6)
                  rv = .true.
                elseif ( i .eq. 7 ) then ! dlw
                  read(buff,*,iostat=ierr) field, XC(nn,7)
                  rv = .true.
                elseif ( i .eq. 8 ) then ! slw
                  read(buff,*,iostat=ierr) field, XC(nn,8)
                  rv = .true.
                elseif ( i .eq. 9 ) then ! bailai
                  read(buff,*,iostat=ierr) field, XC(nn,9)
                  rv = .true.
                elseif ( i .eq. 10 ) then ! tdp 
                  read(buff,*,iostat=ierr) field, XC(nn,10)
                  rv = .true.
                elseif ( i .eq. 11 ) then ! ssc
                  read(buff,*,iostat=ierr) field, XC(nn,11)
                  rv = .true.
                elseif ( i .eq. 12 ) then ! shl
                  read(buff,*,iostat=ierr) field, XC(nn,13)
                  rv = .true.
                elseif ( i .eq. 13 ) then ! br_refl
                  read(buff,*,iostat=ierr) field, rtemp
c                 test if a number is directly after the keyword an undocumented option
                  if (  ierr .ne. 0 ) read(buff,*,iostat=ierr) 
     &              field, st, rtemp
c                   if there was an error, the word 'fixed' is between tr_refl and value 
                  do iw = IC(10), IC(11)
                    SPCIN( iw , 3*nclmax+7+nn ) = rtemp
                  enddo
                  rv = .true.
                elseif ( i .eq. 14 ) then ! tr_refl
c                 test if a number is directly after the keyword: an undocumented option
                  read(buff,*,iostat=ierr) field, rtemp
                  if (  ierr .ne. 0 ) read(buff,*,iostat=ierr) 
     &              field, st, rtemp
c                   if there was an error, the word 'fixed' is between tr_refl and value
                  do iw = IC(10), IC(11)
                    SPCIN( iw , 4*nclmax+7+nn ) = rtemp
                  enddo
                  rv = .true.
                elseif ( i .eq. 15 ) then ! species
                  read(buff,*,iostat=ierr) field, DESC(nn+1)
                  rv = .true.
                elseif ( i .eq. 16 ) then !  leafmodel
c                 if leafmodel = fixed, read reflectance and transmittance
                  call splil2( buff, st, sn )
                  if ( st .eq. 'fixed' ) then
                    read(buff,*,iostat=ierr) field, st, rtemp, ttemp
                    do iw = IC(10), IC(11)
                      SPCIN(iw, 7+nn ) = rtemp
                      SPCIN(iw, nclmax+nn+7) = ttemp
                    enddo
                    rv = .true.
                  endif
                elseif ( i .eq. 17 ) then ! wax_cf
                  read(buff,*,iostat=ierr) field, XC(nn,12)
                  rv = .true.
                elseif ( i .eq. 18 ) then ! use_p
                  read(buff,*,iostat=ierr) field, use_p
                  rv = .true.
                endif
              endif
              if ( rv ) then
c               error checking
                if ( ierr .ne. 0 ) then
                  write(*,*) cln,' Warning: Subroutine rd_ell'//
     &              ' could not read ', field
                  write(*,*) buff 
                  cf_warn = .true.
                else
                  filled(i) = .true.
                endif
              else
c               now fields referring to another section
c               ( i must be 16 ), leafmodel (PROSPECT or LIBERTY)
c               THIS IS CURRENTLY NOT IMPLEMENTED, spectra are read from files (2022)
                if ( sn .ne. ' ' ) then
c                 add to info string calling SECTION type and number of canopy class
                  if ( itype .eq. 1 ) then
                    stemp2 = 'ellipsoid nn='
                    i1 = 13
                  elseif ( itype .eq. 2 ) then
                    stemp2 = 'cone nn='
                    i1 = 8
                  endif
                  write(stemp,'(a,i2,a)') stemp2(1:i1),nn,'  '
                  call addreq( st, sn, stemp )
                  filled(i) = .true.
                else
                  write(*,*) cln,' Warning: Subroutine rd_ell'//
     &              ' could not read ', field
                  cf_warn = .true.
                endif
              endif
c             break out of the "do i=1,nflds" loop
              exit
            endif
          enddo
          if (i .gt. nflds) then
            write(*,*) cln,' WARNING: unknown field: ', field
            cf_warn = .true.
          endif
        endif
      enddo
c     END  ----  do while (loop)
c     check whether there's something missing
      do i = 1, nflds
        if (.not. filled(i) ) then
          if ( .not. ( itype .eq. 1 .and. i .eq. 4 ) ) then
              write(*,*) 'ERROR: required field missing: ',fields(i)
              cf_err = .true.
          endif
        endif
      enddo
c
      if (use_p) then
        p = 1. - XC(nn,11)
        do iw = 1, nspchnl
          if ( iw .lt. IC(10) .or. iw .gt. IC(11) ) then
            SPCIN(iw, 7+nn ) = 0
            SPCIN(iw, nclmax+nn+7) = 0
          else
c             calculate leaf albedo
            w = SPCIN(iw,7+nn)+SPCIN(iw,nclmax+nn+7)
c             preserve R/T ratio. Currently, no better option.
            SPCIN(iw,7+nn) = SPCIN(iw,7+nn)*(1-p)/(1-p*w)
            SPCIN(iw,nclmax+nn+7)=SPCIN(iw,nclmax+nn+7)*(1-p)/(1-p*w)
          endif
        enddo
        write(*,*) 'Scaled leaf R & T using p=',p
      endif
c
      IC(1) = nn
      return
      end ! subroutine rd_ell
c
c
      subroutine rd_frt1( ifh, XC, IC, SPCIN )
c
c     the actual input if SECTION is frt, phase 1: read wavelength range
c  IN: ifh: input file handle,
c  OUT: XC, IC, SPCIN: model input matrices
C   note: XC is currently not used (2022)
      implicit none
      include 'cf_new.h'
      double precision XC(nclmax,*), SPCIN(nspchnl,*)
      integer ifh, IC(*)

      character*260 buff
      integer i, nflds, ierr, ierr2
      parameter (nflds=1)
c       flags
      logical loop, ll, rv 
c       field names
      character*16 fields( nflds ) 
c       whether a value has been found in the file
      logical filled( nflds ) 
c       whether a value is compulsory
      logical reqd( nflds ) 
c       whether to add a line into reqsec
      logical l_add 
      character*16 field
c       requested SECTION type
      character*30 st, sn 
c       max. number of channels ( -1 for any number)
      integer maxchnl 

      common /cfm/ maxchnl
      save /cfm/
c
      data fields /'wavelength'/
      data filled / nflds * .false. /
      data reqd   / nflds * .true. /

      l_add = .true. 
c     -----------------------------------------------------------------
c     LOOP for reading lines from input file
      loop = .true.
      ierr = 0
      do while (loop)
        read (ifh, '(A)', IOSTAT=ierr) buff
        cln = cln + 1
c       error check
        if ( ierr .ne. 0 ) then
          write (*,*) cln,' WARNING: unexpected end of file'
c         exit the loop immediately
c         variables will be checked after the loop, so if all
c         important fields were filled, there is no reason to
c         create an error
          loop = .false.
          exit
        else 
c         strip blanks and check line
          call chsl( ifh, buff, ll, loop )
        endif
        if ( ll ) then
          read(buff,*) field
c         read the value into appropriate global variable
          do i=1,nflds
C           The loop is totally unnecessary, it's here to be in line with other rd_* subroutines          
            rv = .false.
            if ( fields(i) .eq. field ) then
c               determine section type and name
                call splil2( buff, st, sn )
                if ( st .ne. ' ' ) then
                  if ( i .eq. 1 ) then ! wavelength
c                   try to read it to see if it's a single value  
c                   two options: either the value is directly behind  
c                   the word 'wavelength' or the word 'discrete'   
c                   is between them
                    read(st,*,iostat=ierr) SPCIN(1,1)
                    if ( SPCIN(1,1) .gt. 0 ) then
                      IC(10) = 1
                      IC(11) = 1
                      wl_set = .true.
                      l_add = .false.
                    elseif ( st .eq. 'discrete' ) then  
c                     two options: sn contains wavelength in nm or section name
                      read(sn,*,iostat=ierr2) SPCIN(1,1)
                      if ( ierr2 .eq. 0 ) then ! sn contains wavelength 
                        IC(10) = 1
                        IC(11) = 1
                        wl_set = .true.
                        l_add = .false.
                      endif
                    endif
                  endif
c                 a new section has to be read, add it to reqseq
                  if (l_add) call addreq( st, sn, 'frt' )
                  filled(i) = .true.
              endif
c             break out of the "do i=1,nflds" loop
              exit
            endif
          enddo
        endif
      enddo
c     END  ----  do while (loop)
c     check whether there's something missing
      do i = 1,nflds
        if ( reqd(i) .and. .not. filled(i) ) then
          write(*,*) 'ERROR: required field missing: ',fields(i)
          cf_err = .true.
        endif
      enddo
c
c     add the need to re-read frt to get the rest of the data
c     if wl_set .eq. .true., rd_frt2 will be called. If  the required wavelength
c          section is present in the file, it will be found before SECTION FRT as a whole loop
c          in the file is required before FRT will be read again.
      call addreq( 'frt', '', 'frt1' )
c
      return
      end ! subroutine rd_frt1
c
c
      subroutine rd_frt2( ifh, XC, IC, SPCIN, DESC )
c     the actual input if SECTION is frt
c      read everything except the wavelengths used (which are read in rd_frt1)
c  IN: ifh: input file handle,
c  INOUT: IC, SPCIN: model input data matrices
c  OUT: XC: model input data matrices
c    DESC: character array for stand description and class species' names
      implicit none
      include 'cf_new.h'
      double precision XC(nclmax,*), SPCIN(nspchnl,*)
      integer ifh, IC(*)
      character*(*) DESC(*)

      character*260 buff
      integer nflds, ierr, ierr2, iw
      parameter (nflds=14)
c       flags
      logical loop, ll, rv 
c       field names
      character*16 fields( nflds ) 
c       whether a value has been found in the file
      logical filled( nflds ) 
c       whether a value is compulsory
      logical reqd( nflds ) 
c       whether to add a line into reqsec
      logical l_add 
      character*16 field
c       requested SECTION type
      character*30 st, sn 
      double precision tmp_dbl, dage, mindist
      integer i, j, i3, maxchnl

      common /cfm/ maxchnl
      save /cfm/
c
      data fields /'name', 'age', 'refr_idx', 'job_id',
     & 'treeclass', 'groundmodel', 'nlayers', 'skymodel',
     & 'invmodel', 'wavelength', 'angles', 'cubature',
     &  'RAMIbase', 'correction_wl' /
      data filled / nflds * .false. /
      data reqd   / nflds * .true. /
c
c     deal with optional fields: nothing bad happens, if they are not filled:
      reqd(1)=.false.
      reqd(2)=.false.
      reqd(9)=.false.
      reqd(13)=.false.
      reqd(14)=.false.
      dage = -1
c     initialize correction wavelength to beyond index range: no correction
c     this takes effect if the field contains rubbish or is missing
      IC(2) = IC(11)+100
c     -----------------------------------------------------------------
c     LOOP for reading lines from input file
      loop = .true.
      ierr = 0
      do while (loop)
        read (ifh, '(A)', IOSTAT=ierr) buff
        cln = cln + 1
c       error check
        if ( ierr .ne. 0 ) then
          write (*,*) cln,' WARNING: unexpected end of file'
c         exit the loop immediately
c         variables will be checked after the loop, so if all
c         important fields were filled, there is no reason to
c         create an error
          loop = .false.
          exit
        else 
c         strip blanks and check line
          call chsl( ifh, buff, ll, loop )
        endif
        if ( ll ) then
          read(buff,*) field
c         read the value into appropriate global variable
          do i=1,nflds
            rv = .false.
            if ( fields(i) .eq. field ) then
c
c             first requests to read a file
              call splil2( buff, st, sn )
              if ( st .eq. 'file' ) then
                if ( i .eq. 3 ) then ! refr_idx
                  call rfile(sn, SPCIN(1,1), IC(10), IC(11),
     &              SPCIN(1,7), 1, ierr)
                  rv = .true.
                endif
              else
c
c             now data fields
                if ( i .eq. 1 ) then ! avoiding calculated goto
                  read(buff,*,iostat=ierr) field, DESC(1)
                  rv = .true.
                elseif ( i .eq. 2 ) then ! age
                  read(buff,*,iostat=ierr) field, XC(1,25)
                  rv = .true.
                elseif ( i .eq. 3 ) then ! refr_idx
                  read(buff,*,iostat=ierr) field, SPCIN(1,7)
                  do iw = IC(10), IC(11)
                    SPCIN(iw,7) = SPCIN(1,7)
                  enddo
                  rv = .true.
                elseif ( i .eq. 4 ) then ! job_id
                  read(buff,*,iostat=ierr) field, IC(12+nclmax)
                  rv = .true.
                elseif ( i .eq. 7 ) then !  nlayers
                  read(buff,*,iostat=ierr) field, IC(3)
                  rv = .true.
                elseif ( i .eq. 13 ) then ! RAMIbase
                  write(*,*) 'WARNING: deprecated field ',
     &              'RAMIbase in config file, ignoring.'
                  rv = .true.
                elseif ( i .eq. 14 ) then ! correction_wl
                  read(buff,*,iostat=ierr) field, tmp_dbl
                  if ( ierr .eq. 0 ) then
c                     negative values indicate an index in wl array
c                     if index is beyond wl array upper bound -> no correction
                    if ( tmp_dbl .lt. 0 ) then
                      IC(2) = floor( -tmp_dbl )
                    elseif ( tmp_dbl .eq. 0 ) then
c                     zero indicates that correction is computed for each wl separately    
                      IC(2) = 0
                    else
c                     find closest wl, search in all SPCIN(:,1) regardless of lower boundary
                      mindist = 1e8
                      do j=1,IC(11)
                        if (abs(SPCIN(j,1)-tmp_dbl) .lt. mindist) then
                          IC(2)=j
                          mindist=abs(SPCIN(j,1)-tmp_dbl)
                        endif
                      enddo
                    endif
                  endif
                  rv = .true.
                endif
              endif
c             error checking
              if ( rv ) then
                if ( ierr .ne. 0 ) then
                  write(*,*) cln,' Warning: Subroutine rd_frt'//
     &              ' could not read ', field
                  cf_warn = .true.
                else
                  filled(i) = .true.
                endif
              else
c
c               now fields referring to another section
c                ( i must be 5, 6, 8, 9, 10, 11)
c               can't use simple read since the number of fields
c               is unknown ( 1 or 2 )
                l_add = .true.   
c               determine section type and name
                call splil2( buff, st, sn )
                if ( st .ne. ' ' ) then
                  if ( i .eq. 6 ) then ! groundmodel
c                  set ground model type
                    if ( st .eq. 'lambert' ) then
c                     try to read also the reflectance value
                      read(sn,*,iostat=ierr2) tmp_dbl
                      if ( ierr2 .eq. 0 ) then     
c                       constant ground albedo value
                        l_add = .false.
                        do i3 = IC(10), IC(11)
                            SPCIN(i3,2) = tmp_dbl
                            SPCIN(i3,3) = tmp_dbl
                            SPCIN(i3,4) = tmp_dbl
                            SPCIN(i3,5) = tmp_dbl
                        enddo
                      endif
                    else
                      write(*,*) 'ERROR: unknown ground model ',st
                      cf_err = .true.
                    endif
                  elseif ( i .eq. 8 ) then ! skymodel
                    if ( st .eq. 'difsky' ) then
c                     try to read it to/ see if it's a single value, S/Q ratio
                      tmp_dbl = 0
c                     apparently, if read fails, tmp_dbl may remain uninitialized
c                     this is still relevant if 
                      read(sn,*,iostat=ierr) tmp_dbl
c                      for some reason, ierr was not used efore 2022
c                      if ( tmp_dbl .gt. 1e-100 .or.
c     &    ( tmp_dbl .eq. 0 .and. sn(1:1) .eq. '0' ) ) then
                      if (ierr .eq. 0) then
                        do i3 = IC(10),IC(11)
                          SPCIN(i3,6) = tmp_dbl
                        enddo
                        l_add = .false.
                      endif
                    endif
                  elseif ( i .eq. 10 ) then ! wavelength
c                   this was already handled in rd_frt1, do nothing here
c                   this section name is retained to avoid warning of "unknown field"
                    l_add = .false.
                    filled(i) = .true.
                  elseif ( i .eq. 12 ) then ! cubature
c                   try to read it to see if it's 2 numbers
c                   two options: either the numbers are directly
c                   behind the word 'cubature' or the word
c                   'GLcubature' is between 'cubature' & numbers
                    read(st,*,iostat=ierr) IC(13+nclmax)
                    read(sn,*,iostat=ierr)  IC(14+nclmax)
                    if (  IC(13+nclmax) .gt. 0 
     &               .and.  IC(14+nclmax) .gt. 0 ) then
                      l_add = .false.
                    elseif ( st .eq. 'GLcubature' ) then
                      read (buff,*,iostat=ierr) 
     &                  field,st, IC(13+nclmax), IC(14+nclmax)
                      if (  IC(13+nclmax) .gt. 0 
     &                   .and.  IC(14+nclmax) .gt. 0 ) then
                        l_add = .false.
                      endif
                    endif
                  endif
c                 a new section has to be read, add it to reqseq
                  if (l_add) call addreq( st, sn, 'frt' )
                  filled(i) = .true.
                else 
                  write(*,*) cln,' Warning: Subroutine rd_frt'//
     &              ' could not read ', field
                  cf_warn = .true.
                endif
              endif
c             break out of the "do i=1,nflds" loop
              exit
            endif
          enddo
          if (i .gt. nflds) then
            write(*,*) cln,' WARNING: unknown field: ', field
            cf_warn = .true.
          endif
        endif
      enddo
c     END  ----  do while (loop)
c     check whether there's something missing
      do i = 1,nflds
        if ( reqd(i) .and. .not. filled(i) ) then
          write(*,*) 'ERROR: required field missing: ',fields(i)
          cf_err = .true.
        endif
      enddo

      return
      end ! subroutine rd_frt2
c
c
      subroutine isreq( sectyp, secnam, nx )
c     looks the section up in reqsec
c  IN:  sectyp, secnam: section type and name
c  OUT: nx: the index of the section in reqseq (-1 if it is not found)
      implicit none
      include 'cf_new.h'
      character*(*) sectyp, secnam
      integer nx, i  
c
      nx = -1
      do i=1,nsecm
         if ( reqsec(i,1) .eq. sectyp ) then
          if ( (reqsec(i,2)(1:1) .eq. ' ')
     &      .or. (reqsec(i,2) .eq. secnam) ) then
            nx = i
            exit
          endif
        endif
      enddo
      end  ! subroutine isreq
c
c
      subroutine addreq( sectyp, secnam, secpar )
c
c     adds a request for a SECTION to be found in input file
c     into the first empty slot of reqsec
c  IN: the variables forming the row of reqsec
c
      implicit none
      include 'cf_new.h'
      character*(*) sectyp, secnam, secpar
      integer i
c
      do i=1,nsecm
        if ( reqsec(i,1) .eq. ' ' ) then
          reqsec(i,1) = sectyp
          reqsec(i,2) = secnam
          reqsec(i,3) = secpar
          exit
        endif
      enddo
c
      return
      end ! subroutine addreq
c
c
      subroutine remreq( ii )
c
c     removes a request for a SECTION to be found in input file
c     from the slot (=row) ii of reqsec
c  IN: ii :reqsec row number
c
      include 'cf_new.h'
c
      reqsec(ii,1) = ' '
      reqsec(ii,2) = ' '
      reqsec(ii,3) = ' '
c
      return
      end ! subroutine addreq
c
c
      subroutine chsl( ifh, buff, ll, loop )
c
c     checks the line to see if there should be data
c     this is used if we're going through the lines of a section
c   IN: ifh: input file handle
c   IN/OUT: buff: line from input file (will be modified) DECLARED IN cf_new.h
c   OUT: ll: whether line can contain data
c     loop: whether this will NOT end the SECTION
c
      implicit none
      include 'cf_new.h'
      integer ifh
      character*(*) buff
      logical ll, loop
      integer nb
c
      ll = .true.
      loop = .true.
      call fndnb( buff, nb )
      if ( nb .ge. 1 ) then 
c       this is not a blank line
c       check is the section ends or is there more data
        if ( buff(1:4) .eq. 'END ' ) then
          loop = .false.
          ll = .false.
        elseif ( buff(1:7) .eq. 'SECTION ' ) then
          write(*,*) cln,
     &     ' WARNING: new SECTION before previous ended'
          backspace ifh
          cf_warn = .true.
          loop = .false.
          ll = .false.
        elseif ( buff(1:1) .eq. '#' ) then
c          actually, all comments should be stripped before
c            this point
          ll = .false.
        endif 
c     otherwise it has to be data
      else
        ll = .false.
      endif
c
      return
      end ! subroutine chsl
c
c
      subroutine fndlnb( buff , nb )
c
c     gives the location of the last non-blank character in string buff
c  IN:  buff: input string
c  OUT: nb: location of 1st nonblank (-1 if not found)
      implicit none
      character*(*) buff
      integer nb

      nb = len( buff )
 
      do while ( buff(nb:nb) .eq. ' ' ) 
        if ( nb .gt. 1 ) then
          nb = nb - 1
        else
          nb = -1
          exit
        endif
      enddo
      return
      end ! subroutine fndlnb
c
c 
      subroutine rfile( fname, X, n1, n2, parm, n, ierr )
c
c     fill matrix parm from file fname
c     also can modify the wavelength limits of the model
c   IN: fname: file name
c     X : array of wavelength values
c     n: no. of columns in parm
c   INOUT:  n1, n2: indices in X giving the used wavelength range
c        commonly, these values are stored in IC(10) and IC(11), respectively
c   OUT: parm: parameter matix
c     ierr: error flag
      implicit none
      include 'cf_new.h'
      character*(*) fname
      integer n1, n2, n, ierr, ierr2
      double precision X(nspchnl)
      double precision parm(nspchnl,n)
c
      integer i, j, jl, jp, ndcount, nb, ire
      integer bestwl
      double precision bestdist, tempdist
      character*260 buffx
c     buffers for wavelengths and spectra as read from file
c        note: these can be larger than what's allowed in the model, so 
c        only a subset is used in the model. Allow ample space for storage
      double precision B_X(10*nspchnl)
      double precision B_parm(10*nspchnl,n)
      logical n12changed
      logical datafound
c     datafound indicates if USEFUL data has been found, corresponding to 
c       a wavelength in X
c
      ire = 99
      ierr = 0
      jl = 0 ! line number in file
      jp = 1 ! location in X and parm
      ndcount = 0 ! non-data line count
c
      open (ire,file=fname,status='old',iostat=ierr)
      if ( ierr .gt. 0 ) then
        write(*,*) 'ERROR opening file ',fname
      else      
        rewind ire
c       LOOP OVER FILE and read its contents into B_X and B_parm
 110    continue   
c         read a line into buffer
          read (ire, '(A)', iostat=ierr) buffx
c         when reading the last line, read returns err=-1
c         regardless of whether this last line contains any data or not
          jl = jl + 1
c        remove comments and leading blanks
          call fndnb( buffx , nb )
          if ( nb .gt. 0 ) then
            read(buffx,*,iostat=ierr2) B_X(jp), ( B_parm(jp,j), j=1,n )
            if ( ierr2 .gt. 0 ) then
              ndcount = ndcount + 1 
            else
              jp = jp + 1
            endif
          endif
        if ( ierr .eq. 0 .and. jp .lt. 10*nspchnl ) goto 110
c       ---------------- end of until-loop 
c       rewind jp to last added data location
        jp = jp - 1
        write(*,*) "rfile: scanned ",jp," and skipped", ndcount,
     &  " lines in ",trim(fname)
        if ( jp .gt. 0 ) ierr = 0

c       match model spectra values in X and file spectra values in B_X
c       a good old double loop
        n12changed = .false.
        datafound = .false.
        do i=n1,n2
          bestwl = -1
          bestdist = 10*wldelta  
          do j=1,jp
            tempdist = abs(X(i)-B_X(j))
            if ( tempdist .lt. bestdist ) then
                bestdist = tempdist
                bestwl = j
            endif
          enddo

          if ( bestdist .lt. wldelta ) then
            datafound = .true.
            do j=1,n
              parm(i,j)=B_parm(bestwl,j)
            enddo
          else
c           no data at this wavelength, 
c               increase spectrum start wavelength if no data found, or 
c               decrease spectrum end wavelength if we have some data already
            n12changed = .true.
            if ( datafound ) then
              n2=i-1
c             spectrum needs to be continuous, skip further wavelengths
              exit
            else
              n1=i+1
            endif
          endif
        enddo
c       check if we have changed wl limits of the model
        if (n12changed) then
          if (n1 .gt. n2) then
            write(*,*) "STOPPING: no suitable spectral data found in ",
     &       fname
            ierr=100
          else
            write(*,"(A,F6.1,A,F6.1,2A)") "changed spectrum range to ",
     &        X(n1),"--",X(n2)," due to limits in file ", trim(fname)
          endif
        endif
c         set the error flag to 0 if loop was finished due to EOF
c           or at least some values were read (read sets it to -1)
  
c
      endif
      close (ire)
c
      return
      end  ! subroutine rfile
