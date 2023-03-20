*
*     cf_new.h header file for functions reading the new input format
*
*
*     defines common /cf_new/
*     and also some generally used non-global variables
*
      include 'frtpar.h'
      
      integer nsecm
      parameter ( nsecm = 15 ) ! max. number of sections in input file
      character*32 reqsec(nsecm,3) ! which sections are required
!        columns: type | name | the name of the SECTION requiring it
      double precision wldelta 
      logical cf_warn ! whether to report of a warning
      logical cf_err ! whether to report an error
      logical wl_set ! whether wavelength range has been read into SPCIN(:,1)
      integer cln ! current line number in input file
      integer nn ! the current tree class 
      common /cf_new/ reqsec, wldelta,
     &   cln,nn, cf_warn, cf_err, wl_set
      save /cf_new/
*
*  nsecm:  max. number of sections in input file
*  cminwl: minimum wavelength model can use
*  cmaxwl: max. wavelength model can use
*  wldelta: wavelength accuracy for matching different spectral files: the largest tolerated wl difference [nm]
*  min/max wavelength are determined from input files that contain e.g., measured spectra
