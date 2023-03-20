# functions used by FRT to compute reflectance and transmittance
#      versions which use modules written in f77
#      either this or the corresponding frtfunctions_py  needs to be imported by frt
#  f77 modules need to be compiled with f2py and require a fortran 90 compiler (e.g., gortran)
#  sample commands for windows are given within each modue, but they depend on the specific fortran installation, e.g.,
#    f2py.exe -c --compiler=mingw32 -m pi11d pi11d.f
#  alternatively, fortran can be compiled from within python:
#  numpy.f2py.compile(source, modulename='untitled', extra_args='', verbose=True, source_fn=None, extension='.f')

import numpy as np
from frtfunctions import *
# frt fortran77 modules:
import spooi
import enel3
import bck3
# fill in the common blocks in the fortran77 modules
spooi.pidr.pi = np.pi
spooi.pidr.dr = np.pi/180
enel3.pidr.pi = np.pi
enel3.pidr.dr = np.pi/180
bck3.pidr.pi = np.pi
bck3.pidr.dr = np.pi/180

# NOTE: enel3 also uses the /volint/ common block, filled in elsewhere
# load the /frtpar/ common block
enel3.enel_incl()
# /frtpar variables can now be accessed as e.g. enel3.frtpar.fnclmax[()]

def hetk8s_singledirection(theta_vec, l_elli, shl, stdns, htr, hc1, hc2, rcr, dbh, glmp, ulg):
    """ Calculation of unidirectional upward gap probabilities in crown layer

    Tõenäosuste arvutamine               A. Kuusk   1.12.2000
        Modified by Matti Mõttus (2020). Basically, a wrapper around spooi.
        spooi is (currently?) implemented in fortran77 and imported as library

    hetk8sA() in fortran code with reordered input parameters
    Args: except theta_vec, all arrays of the same length (one value for tree class)
    theta_vec: np.ndarray of zenith angles
    l_elli: whether crown shape is ellipsoid (logical)
    shl: shoot length
    stdns: stand density [1/m]
    htr: tree height [m]
    hc1: crown length, ellipse / conical part [m]
    hc2: crown length, cylinder [m]
    rcr: crown radius [m]
    dbh: trunk diameter at brest height [m]
    glmp: Fischer's grouping index
    ulg: uuu / 2 (=uuu*G?)???, uuu: effective plant area density (leaves corrected for clumping + branches)

    Returns:
    gaps: gap probability (np.ndarray)
    """

    ncl = len(l_elli )

    gaps = np.zeros_like(theta_vec)

    # these should be the coordinates of rays entering and exiting the crown
    #  used only for hotspot computations when integrting over tree crowns
    # (A. Kuusk, The hot spot effect in plant canopy reflectance. In R. Myneni
    # and J. Ross (Eds.), Photon-Vegetation Interactions. Springer, Berlin, 1991, 139-159.
    # set here to zero

    x1=y1=z1=x2=y2=z2=x3=y3=z3=0.0

    #                                  gap probability on ground
    #                                        gaps(maapinnal)
    # additional hotspot parameters
    rlls1 = rllv1 = rllv3 = chs3 = sphi2 = 0.0
    # angles between view and sun directions -- not required for unidirectional computations
    cphi2  = 1
    chs1   = 0
    calph = 1
    # f2py requires input arrays to be of exactly correct shape
    nclmax = enel3.frtpar.fnclmax[()]
    l_elli_F = np.zeros(nclmax, dtype=bool)
    shl_F = np.zeros(nclmax)
    stdns_F = np.zeros(nclmax)
    htr_F = np.zeros(nclmax)
    dbh_F = np.zeros(nclmax)
    hc1_F = np.zeros(nclmax)
    hc2_F = np.zeros(nclmax)
    rcr_F = np.zeros(nclmax)
    glmp_F = np.zeros(nclmax)
    ulg_F = np.zeros(nclmax)
    ncl = len(l_elli)
    l_elli_F[0:ncl] = l_elli
    shl_F[0:ncl] = shl
    stdns_F[0:ncl] = stdns
    htr_F[0:ncl] = htr
    dbh_F[0:ncl] = dbh
    hc1_F[0:ncl] = hc1
    hc2_F[0:ncl] = hc2
    rcr_F[0:ncl] = rcr
    glmp_F[0:ncl] = glmp
    ulg_F[0:ncl] = ulg

    for i,t in enumerate(theta_vec):
        stheta = np.sin(t)
        ctheta = np.cos(t)

        #   calculate aas: gap probability, Sun direction
        aasi, poodi = spooi.spooi( l_elli_F, ncl, ulg_F, shl_F, stdns_F, htr_F,
            hc1_F, hc2_F, rcr_F, dbh_F, glmp_F,
            x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
            stheta, ctheta, stheta, ctheta, sphi2, cphi2, calph, chs1, chs3)
        gaps[i] = aasi

    return gaps  #  --------------------  hetk8s_singledirection

def hetk8s_bidirectional(thetv_vec, phi_vec, thets, nzs,
    l_elli, shl, stdns, htr, hc1, hc2, rcr, dbh, glmp, ulg):
    """ Calculation of bidirectional (view, sun) gap probabilities in crown layer.

    Based on Fortran77 code Tõenäosuste arvutamine, A. Kuusk 1.12.2000
    Modified by Matti Mõttus (2020, 2022) to include also a loop over azimuth
    part of hetk8sB
     all directions ( thetv_vec x phi_vec ) are sampled

    Args:
    thetv_vec: vector of view zenith angles [rad], e.g. knots of G-L quadrature (zenith angle)
    phi_vec: vector of relative view azimuth angles for the quadrature [rad]
    thets: sun zenith angle [rad]
    l_elli (logical): whether crown shape is ellipsoid (vector of length ncl)
    shl: shoot length (vector of length ncl)
    stdns: stand density [1/m] (vector of length ncl)
    htr: tree height [m] (vector of length ncl)
    hc1: crown length, ell / con [m] (vector of length ncl)
    hc2: crown length, cylinder [m] (vector of length ncl)
    rcr: crown radius [m] (vector of length ncl)
    dbh: trunk diameter at brest height [m] (vector of length ncl)
    glmp: Fischer's grouping index (vector of length ncl)
    ulg: uuu / 2 (=uuu*G?)???,  uuu: effective plant area density (leaves corrected for clumping + branches) (vector of length ncl)

    Returns:
    psgvuQ: bidirectional (sun, view) gap probability, upper hemisphere. 2-dim array for the quadrature
    """
    # f2py requires input arrays to be of exactly correct shape
    nclmax = enel3.frtpar.fnclmax[()]
    l_elli_F = np.zeros(nclmax, dtype=bool)
    shl_F = np.zeros(nclmax)
    stdns_F = np.zeros(nclmax)
    htr_F = np.zeros(nclmax)
    hc1_F = np.zeros(nclmax)
    hc2_F = np.zeros(nclmax)
    rcr_F = np.zeros(nclmax)
    dbh_F = np.zeros(nclmax)
    glmp_F = np.zeros(nclmax)
    ulg_F = np.zeros(nclmax)
    ncl = len(l_elli)
    l_elli_F[0:ncl] = l_elli
    shl_F[0:ncl] = shl
    stdns_F[0:ncl] = stdns
    htr_F[0:ncl] = htr
    hc1_F[0:ncl] = hc1
    hc2_F[0:ncl] = hc2
    rcr_F[0:ncl] = rcr
    dbh_F[0:ncl] = dbh
    glmp_F[0:ncl] = glmp
    ulg_F[0:ncl] = ulg

    # output matrix
    psgvuQ = np.zeros( ( len(thetv_vec), len(phi_vec) ) )

    x1 = y1 = z1 = x2 = y2 = z2 = x3 = y3 = z3 = 0
    # rlls1, rllv1, rllv3, chs1 and chs3 are between-crown hotspot parameters in spooi.f
    rlls1 = rllv1 = rllv3 = chs1 = chs3 = 0
    stheta1 = np.sin(thets)
    ctheta1 = np.cos(thets)
    for iph,phi in enumerate( phi_vec ):
        sphi   = np.sin(phi)
        cphi   = np.cos(phi)
        for ith,theta2 in enumerate(thetv_vec):
            # p(sunlit, ground) and pooui(ground)
            # p(sunlit, ground) ja pooui(maapinnal)
            stheta2 = np.sin(theta2)
            ctheta2 = np.cos(theta2)
            calph  = stheta2*stheta1*cphi + ctheta2*ctheta1
            # calph=cos(alpha), alpha is the angle between Sun & view directions

            pooui, poodi = spooi.spooi(l_elli_F, ncl, ulg_F, shl_F, stdns_F,
                htr_F, hc1_F, hc2_F, rcr_F, dbh_F, glmp_F,
                x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
                stheta1, ctheta1, stheta2, ctheta2, sphi, cphi, calph, chs1, chs3)

            psgvuQ[ith,iph] = pooui
        #
    return psgvuQ   # --------------------  hetk8s_bidirectional

def hetk8s_integrated(thetv_vec, phi_vec, thets, nzs,
    l_elli, shl, stdns, htr, hc1, hc2, rcr, dbh, glmp, ulg, uuu):
    """ Calculation of the bidirectional gap probabilities integrated over the crowns.

    Based on Fortran77 code Tõenäosuste arvutamine, A. Kuusk 1.12.2000
    Modified by Matti Mõttus (2022) to include loops over both zenith and azimuth
    A wrapper around the original f77 enel imported as a module
    part of hetk8sB
     all directions ( thetv_vec x phi_vec ) are sampled

    Args:
    thetv_vec: vector of view zenith angles [rad], e.g. knots of G-L quadrature (zenith angle)
    phi_vec: vector of relative view azimuth angles for the quadrature [rad]
    thets: sun zenith angle [rad]
    nzs: number of crown layers in numerical integration
    l_elli (logical): whether crown shape is ellipsoid (vector of length ncl)
    shl: shoot length (vector of length ncl)
    stdns: stand density [1/m] (vector of length ncl)
    htr: tree height [m] (vector of length ncl)
    hc1: crown length, ell / con [m] (vector of length ncl)
    hc2: crown length, cylinder [m] (vector of length ncl)
    rcr: crown radius [m] (vector of length ncl)
    dbh: trunk diameter at brest height [m] (vector of length ncl)
    glmp: Fischer's grouping index (vector of length ncl)
    ulg: uuu / 2 (=uuu*G?)??? (vector of length ncl)
    uuu: effective plant area density (leaves corrected for clumping + branches) (vector of length ncl)

    Returns:
    bdgfuT: Integral of bidirectional gap probability over a crown, upward direction (3-dim np.ndarray)
    bdgfdT: Integral of bidirectional gap probability over a crown, downward direction  (3-dim np.ndarray)
    btr1ukT: probability of seeing sunlit trunk from above  (3-dim np.ndarray)
    btr1dkT: probability of seeing sunlit trunk from below  (3-dim np.ndarray)
    Note: compared with original frt arrays, the output arrays were transposed in 2020    and the letter T was added to their names
    ===WARNING== the contents of bdgfuT, bdgfdT, btr1ukT, btr1dkT was guessed from
    comment lines without going through the code line by line!
    """
    # f2py requires input arrays to be of exactly correct shape
    nclmax = enel3.frtpar.fnclmax[()]
    l_elli_F = np.zeros(nclmax, dtype=bool)
    shl_F = np.zeros(nclmax)
    stdns_F = np.zeros(nclmax)
    htr_F = np.zeros(nclmax)
    hc1_F = np.zeros(nclmax)
    hc2_F = np.zeros(nclmax)
    rcr_F = np.zeros(nclmax)
    dbh_F = np.zeros(nclmax)
    glmp_F = np.zeros(nclmax)
    ulg_F = np.zeros(nclmax)
    ncl = len(l_elli)
    l_elli_F[0:ncl] = l_elli
    shl_F[0:ncl] = shl
    stdns_F[0:ncl] = stdns
    htr_F[0:ncl] = htr
    dbh_F[0:ncl] = dbh
    hc1_F[0:ncl] = hc1
    hc2_F[0:ncl] = hc2
    rcr_F[0:ncl] = rcr
    glmp_F[0:ncl] = glmp
    ulg_F[0:ncl] = ulg

    # Fill the knots and weights of cubature on a sphere (for integrating
    # gap probability over a crown ellipsoid) -- used in enel3
    netst, xetst, yetst, zetst, aetst = cubell9_py()
    # make this available to the enel3 fortran77 module
    #   a bit more space is reerved in frtpar.h than actually used
    enel3.volint.netst = netst
    enel3.volint.xetst[0:len(xetst)]  = xetst
    enel3.volint.yetst[0:len(yetst)]  = yetst
    enel3.volint.zetst[0:len(zetst)]  = zetst
    enel3.volint.aetst[0:len(aetst)]  = aetst
    # Knots and weights of cubature on a circle (for integrating
    # gap probability over a cylindrical crown) in the enel3 f77 module
    nctst, xctst, yctst, actst = cubcirc_py()
    enel3.volint.nctst = nctst
    enel3.volint.xctst[0:len(xctst)] = xctst
    enel3.volint.yctst[0:len(yctst)] = yctst
    enel3.volint.actst[0:len(actst)] = actst
    # The points of the quadrature for the vertical direction
    #   between -1 and 1, also for the enel3 f77 module
    zctst, acztst = gauleg_numpy(-1, 1, 2*nzs)
    enel3.volint.zctst[0:len(zctst)] = zctst
    enel3.volint.acztst[0:len(acztst)] = acztst

    nth = len(thetv_vec)
    nph = len(phi_vec)
    bdgfuQ = np.zeros( (ncl,nth,nph) )
    bdgfdQ = np.zeros( (ncl,nth,nph) )
    btr1ukQ = np.zeros( (ncl,nth,nph) )
    btr1dkQ = np.zeros( (ncl,nth,nph) )

    nknotm = enel3.frtpar.fnknotm
    bdgfuT = np.zeros( (nclmax,nknotm), order='F' )
    bdgfdT = np.zeros( (nclmax,nknotm), order='F' )
    btr1ukT = np.zeros( (nclmax,nknotm), order='F' )
    btr1dkT = np.zeros( (nclmax,nknotm), order='F' )

    # start the actual computation (call fortran 77 subroutine with extra loop for azimuth)
    for iph,phi in enumerate(phi_vec):
        for ith,theta2 in enumerate(thetv_vec):
            # subroutine enel in f77 fills only one element in bdgfuT, bdgfdT, btr1ukT and btr1dkT at a time
            for icl in range(ncl):
                ul = uuu[icl]
                # for some strange reason, the code wants uuu as a scalar (the final user is bck3.f)
                #  enel3.enel fills
                enel3.enel( l_elli_F, icl+1, ith+1, ncl, ul, shl_F, stdns_F, htr_F, dbh_F,
                    hc1_F, hc2_F, rcr_F, ulg_F, glmp_F, thets, theta2, phi, nzs,
                    bdgfuT, bdgfdT, btr1ukT, btr1dkT)

        # copy to actual tree class and angle sizes
        bdgfuQ[0:ncl,0:nth,iph] = bdgfuT[0:ncl,0:nth]
        bdgfdQ[0:ncl,0:nth,iph] = bdgfdT[0:ncl,0:nth]
        btr1ukQ[0:ncl,0:nth,iph] = btr1ukT[0:ncl,0:nth]
        btr1dkQ[0:ncl,0:nth,iph] = btr1dkT[0:ncl,0:nth]
    return bdgfuQ, bdgfdQ, btr1ukQ, btr1dkQ  # ----- hetk8s_integrated

