*  frtpar.h
c  include-file for program FRT    A. Kuusk  21.02.2002
*   modified by M. Mottus 2020
*  this file should be include'd before any other header files by all subroutines
 
* First, fixed integers defining maximum array sizes
      integer nclmax, ncub, nspchnl, nknotm, nphim
      parameter (nclmax = 10)
      parameter (ncub = 48)
c NOTE: ncub has to be probably fixed to 48, it is used to dimension /volint/ variables
      parameter (nspchnl = 2001)
      parameter (nknotm = 40 )
      parameter (nphim = 40)

* what was previously in volint.h
      integer netst, nctst
      double precision xetst(ncub), yetst(ncub), zetst(ncub), 
     & aetst(ncub), xctst(ncub), yctst(ncub), zctst(ncub), 
     & actst(ncub), acztst(ncub)
      logical l_volint
      common /volint/ netst, nctst,
     & xetst, yetst, zetst, aetst, xctst, yctst, zctst, actst, acztst,
     & l_volint

      double precision pi, dr
      common /pidr/ pi, dr

*  dr: degrees to radians conversion factor
*  pi: 3.14
*  ncub: no. of quadrature knots. Used in cubell9 & cubcirc11 for integrating over ellipsoids
*    should probably be at least 48, DO NOT CHANGE TO ANYTHING ELSE HERE!
*  nknotm: max. no. of quadrature (over zenith angle) knots
*  nphim: max. no. of quadrature (over azimuth angle) knots
*  nspchnl: number of spectral channels 
*     common /volint/: variables for volume integral over a crown
*  netst
*  nctst
*  xetst
*  yetst
*  zetst
*  aetst
*  xctst
*  yctst
*  zctst
*  actst
*  acztst
*  l_volint: flag whether the values in /volint/ are set
