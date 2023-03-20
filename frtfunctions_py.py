# functions used by FRT to compute reflectance and transmittance
#  pure python implementation
#      either this or the corresponding frtfunctions_f77 needs to be imported by frt
import numpy as np
from frtfunctions import *

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

    for i,t in enumerate(theta_vec):
        stheta = np.sin(t)
        ctheta = np.cos(t)

        #   calculate aas: gap probability, Sun direction
        aasi, poodi = spooi( l_elli, ulg, shl, stdns, htr,
            hc1, hc2, rcr, dbh, glmp,
            x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
            stheta, ctheta, stheta, ctheta, sphi2, cphi2, calph, chs1, chs3)
        gaps[i] = aasi

    return gaps
#   --------------------  end hetk8s_singledirection


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

            pooui, poodi = spooi(l_elli, ulg, shl, stdns,
                htr, hc1, hc2, rcr, dbh, glmp,
                x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
                stheta1, ctheta1, stheta2, ctheta2, sphi, cphi, calph, chs1, chs3)

            psgvuQ[ith,iph] = pooui
        #
    return psgvuQ
#   -------------------- end hetk8s_bidirectional


def hetk8s_integrated(thetv_vec, phi_vec, thets, nzs,
    l_elli, shl, stdns, htr, hc1, hc2, rcr, dbh, glmp, ulg, uuu):
    """ Calculation of the bidirectional gap probabilities integrated over the crowns.

    Based on Fortran77 subroutine enel (from  enel3.f) modified to include the
    loop in hetk8s [vo~ralt hajuv energia     A. Kuusk     29.08.2002]
    Modified by Matti Mõttus (2022) to include loops over both zenith and azimuth
    part of hetk8sB
      thetv_vec and phi_vec are vectors, all directions ( thetv_vec x phi_vec ) are sampled

    Integreerib kahe suuna avatuse tõenäosuse kroonil bgp*pooi,
    3D integraal ellipsoidis arvutatakse Hammer-Stroudi kubatuuriga,
    Hammer, P.C., Stroud, A.H. Numerical evaluation of multiple
    integrals II. Math. Tables and other Aids Comp. 12 (1958), No 64, 272-280
    ja koonuses ning silindris Gauss-Legendre'i kvadratuuriga z-koordinaadi
    järgi ning Mysovskikh (1981) algoritmiga 16.33 x-y-tasandis.
    Kvadratuuride sõlmed ja kaalud on include-failis volint.h.
    Tõenäosused bgp ja pooi arvutatakse bck3-s, integraalid
    salvestatakse massiivi bdgf. Valgustatud tüve nägemise
    tõenäosused salvestatakse massiivi btr1ik.
    Note: compared with original frt arrays, the output arrays were transposed in 2020
       and the letter T was added to their names

    Args:
    thetv_vec: vector of view zenith angles [rad], e.g. knots of G-L quadrature (zenith angle)
    phi_vec: vector of relative view azimuth angles for the quadrature [rad]    thets: sun zenith angle [rad]
    nzs: number of crown layers in numerical integration
    l_elli (logical): whether crown shape is ellipsoid (vector of length ncl)
    shl: shoot length (vector of length ncl)
    stdns: stand density [1/m] (vector of length ncl)
    htr: tree height [m] (vector of length ncl)
    hc1: crown length, ellipsid or cone [m] (vector of length ncl)
    hc2: crown length, cylinder [m] (vector of length ncl)
    rcr: crown radius [m] (vector of length ncl)
    dbh: trunk diameter at brest height [m] (vector of length ncl)
    glmp: Fischer's grouping index (vector of length ncl)
    ulg: uuu / 2 (=uuu*G?)??? (vector of length ncl)
    uuu: effective plant area density (leaves corrected for clumping + branches) (vector of length ncl)

    Returns:
    bdgfuQ: Integral of bidirectional gap probability over a crown, upward direction  (3-dim np.ndarray)
    bdgfdQ: Integral of bidirectional gap probability over a crown, downward direction  (3-dim np.ndarray)
    btr1ukQ: probability of seeing sunlit trunk from above  (3-dim np.ndarray)
    btr1dkQ: probability of seeing sunlit trunk from below  (3-dim np.ndarray)
    """
    ncl = len(l_elli)
    sthets = np.sin(thets)
    cthets = np.cos(thets)

    # Fill the knots and weights of cubature on a sphere (for integrating
    # gap probability over a crown ellipsoid)
    netst, xetst, yetst, zetst, aetst = cubell9_py()
    nctst, xctst, yctst, actst = cubcirc_py()
    zctst, acztst = gauleg_numpy(-1, 1, 2*nzs)

    # output matrices
    nth = len(thetv_vec)
    nph = len(phi_vec)
    bdgfuQ = np.zeros( (ncl,nth,nph) )
    bdgfdQ = np.zeros( (ncl,nth,nph) )
    btr1ukQ = np.zeros( (ncl,nth,nph) )
    btr1dkQ = np.zeros( (ncl,nth,nph) )

    # start the actual computation (translated from fortran 77) with extra loop for azimuth
    for iph,phi in enumerate(phi_vec):
        sphi   = np.sin(phi)
        cphi   = np.cos(phi)
        for ith,thetv in enumerate(thetv_vec):
            sthetv = np.sin(thetv)
            cthetv = np.cos(thetv)
            calph  = sthetv*sthets*cphi + cthetv*cthets
            #alph = np.arccos(calph) on nurk kiirte r1 ja r2 vahel (rv ja rs vahel)
            for icl in range(ncl):
                xi = x2 = x3 = yj = y2 = y3 = rlls1 = rllv1 = rllv3 = chs1 = chs3 = cell = 0.0
                ul = uuu[icl]
                lellips = l_elli[icl]
                aell    = rcr[icl]
                a2      = aell**2
                sumu    = 0.0
                sumd    = 0.0
                htree   = htr[icl]
                ntstz   = 2*nzs
                if lellips: #                                 ellipsoid
                    hkroon = hc1[icl]
                    cell   = hkroon*0.5
                    for itst in range(netst):
                        xtsti = xetst[itst]*aell
                        ytsti = yetst[itst]*aell
                        ztsti = zetst[itst]*cell
                        atsti = aetst[itst]*a2*cell
                        sxui, sxdi = bck3(l_elli, icl, ul, shl, stdns, htr, dbh,
                            hc1, hc2, rcr, ulg, glmp, sthets, cthets, sthetv, cthetv,
                            sphi, cphi, calph, xtsti, ytsti, ztsti)
                        sumu = sumu + sxui*atsti
                        sumd = sumd + sxdi*atsti
                else:   #            not ellipsoid
                    vhelko = hc1[icl]
                    vhcil  = hc2[icl]
                    hkroon = vhelko + vhcil
                    for itst in range(ntstz): #               cylinder
                        if vhcil > 0.0:
                            sxyu   = 0.0
                            sxyd   = 0.0
                            aztsti = acztst[itst]/2*vhcil
                            ztsti  = (zctst[itst] - 1)/2*vhcil
                            rcirc  = aell
                            for jtst in range(nctst):
                                xtsti = xctst[jtst]*rcirc
                                ytsti = yctst[jtst]*rcirc
                                atsti = actst[jtst]*rcirc**2
                                sxui, sxdi = bck3(l_elli, icl, ul, shl, stdns, htr, dbh,
                                    hc1, hc2, rcr, ulg, glmp, sthets, cthets, sthetv, cthetv,
                                    sphi, cphi, calph, xtsti, ytsti, ztsti)
                                sxyu = sxyu + sxui*atsti
                                sxyd = sxyd + sxdi*atsti
                            sumu = sumu + sxyu*aztsti
                            sumd = sumd + sxyd*aztsti

                        if vhelko > 0: #                      cone
                            sxyu   = 0.0
                            sxyd   = 0.0
                            aztsti = acztst[itst]/2*vhelko
                            ztsti  = (zctst[itst] + 1)/2*vhelko
                            rcirc  = aell*(1 - ztsti/vhelko)
                            for jtst in range(nctst):
                                xtsti = xctst[jtst]*rcirc
                                ytsti = yctst[jtst]*rcirc
                                atsti = actst[jtst]*rcirc**2
                                sxui, sxdi = bck3(l_elli, icl, ul, shl, stdns, htr, dbh,
                                    hc1, hc2, rcr, ulg, glmp, sthets, cthets, sthetv, cthetv,
                                    sphi, cphi, calph, xtsti, ytsti, ztsti)
                                sxyu = sxyu + sxui*atsti
                                sxyd = sxyd + sxdi*atsti
                            sumu = sumu + sxyu*aztsti
                            sumd = sumd + sxyd*aztsti
                        #
                    #
                bdgfuQ[icl,ith,iph] = sumu
                bdgfdQ[icl,ith,iph] = sumd
                #
                #    Valgustatud tu"ve na"gemise to~ena"osus / probability of seeing sunlit trunk
                #
                btr1u = 0.0
                btr1d = 0.0
                utyu  = 0.0
                utyd  = 0.0

                tgths = sthets/cthets
                if lellips:
                    zv = np.sqrt(a2 + (cell*tgths)**2)/tgths
                else:
                    zv = aell/tgths
                zvari = htree - hkroon - zv # zvari - krooni varju ko~rgus, height of crown shadow

                if zvari > 0:
                    if thetv < thets:
                        tgthv = sthetv/cthetv
                        if tgthv == 0:
                            zview = 0
                        else:
                            if lellips:
                                zv = np.sqrt(a2 + (cell*tgthv)**2)/tgthv
                            else:
                                zv = aell/tgthv
                            zview = max( (htree - hkroon - zv), 0 )
                    else:
                        zview = zvari
                #               na"htava valgustatud tüve pindala uty
                    utyd = dbh[icl]*np.pi/2.0*zvari
                    utyu = dbh[icl]*np.pi/2.0*zview
                    #                      pooi integreerimine [0, zvari]
                    #                       integration of pooi over [0, zvari]
                    dzenel = hkroon/ntstz  # vertical integration step
                    if zvari <= dzenel:
                        dztrunk = zvari
                    else:
                        dztrunk = dzenel
                    nz = int(zvari/dztrunk + 0.5) # number of integrating steps
                    nview = int(zview/dztrunk + 0.5) + 1
                    dztrunk = zvari/nz # make sure we have the correct value for integration step
                    xi      = 0.0
                    x2      = xi
                    x3      = xi
                    yj      = 0.0
                    y2      = yj
                    y3      = yj
                    rlls1   = 0.0
                    rllv1   = 0.0
                    rllv3   = 0.0
                    chs1    = 0.0
                    chs3    = 0.0

                    if nview < nz: #integrate btr1u only until nview
                        for k in range(1,nz+2):
                            zk   = (k-1)*dztrunk
                            z2   = zk
                            z3   = zk
                            pooui, poodi = spooi(l_elli, ulg, shl, stdns,
                                htr, hc1, hc2, rcr, dbh, glmp,
                                xi, yj, zk, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
                                sthets, cthets, sthetv, cthetv, sphi, cphi, calph, chs1, chs3)
                            if (k==1) or (k==nz+1):
                                    poodi = poodi/2
                            btr1d = btr1d + poodi
                            if (k==1) or (k==nview):
                                    pooui = pooui/2
                            if k <= nview:
                                btr1u = btr1u + pooui
                        btr1d = btr1d/nz
                        btr1u = btr1u/nview
                    else:
                        for k in range(1,nz+2):
                            zk   = (k - 1)*dztrunk
                            z2   = zk
                            z3   = zk
                            pooui, poodi = spooi(l_elli, ulg, shl, stdns,
                                htr, hc1, hc2, rcr, dbh, glmp,
                                xi, yj, zk, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
                                sthets, cthets, sthetv, cthetv, sphi, cphi, calph, chs1, chs3)
                            if (k==1) or (k==nz+1):
                                poodi = poodi/2
                                pooui = pooui/2
                            btr1d = btr1d + poodi
                            btr1u = btr1u + pooui
                        btr1d = btr1d/nz
                        btr1u = btr1u/nz
                    btr1ukQ[icl,ith,iph] = btr1u*utyu
                    btr1dkQ[icl,ith,iph] = btr1d*utyd
                # if zvari > 0
            # close the loop over ncl
        # close the loop over thetv_vec
    # close the loop over phi_vec
    return bdgfuQ, bdgfdQ, btr1ukQ, btr1dkQ
#   --------------------- end hetk8s_integrated


def spooi( lelli, ulg, shl, stdns, htr, hc1, hc2, rcr, dbh,
    glmp, x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
    sthets, cthets, sthetv, cthetv, sphi, cphi, calph, chs1, chs3 ):
    """ bidirectional probability of between-crown gaps
         (calculated using Eq. (5) in manual ver. 05.2005)

      subroutine spooi
     & (lelli, ncl, ulg, shl, stdns, htr, hc1, hc2, rcr, dbh,
     & glmp, x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls1, rllv1, rllv3,
     & sthets, cthets, sthetv, cthetv, sphi, cphi, calph, chs1, chs3,
     & pooui, poodi)

    Kahe suuna avatuse to~ena"osus väljaspool võra.   A. Kuusk 22.05.1998
    pooui - ülemine, poodi - alumine poolsfäär,  A. Kuusk 22.11.2000
    Translation to python: Matti Mõttus 2022

    Kroonide projektsioonid asendatakse sama pindalaga ringidega.
     (xi, yi, zi) on kiirte päikese ning vaatleja suunas üles ja alla
    võrast väljumise punktid: i = 1 - vaatleja suunas üles,
    i = 2  - päike, i = 3 - vaatleja suunas alla

    scomm oli väikese sl8 korral vale.                    25.04.2002

    ground: z=0

    Args:
    Sun and view directions given by the pairs [sin(t),cos(t)]:
      [sthets, cthets] & [sthetv, cthetv]
    if  [sthets, cthets] = [sthetv, cthetv], calculate monodirectional
      transmittance (one direction only, no hotspot, etc.)
    calph=cos(alpha), alpha is the angle between Sun & view directions

    Returns:
    pooui: bidirectional gap probability, upper hemisphere
    poodi: bidirectional gap probability, lower hemisphere
    """
    eps1 = 0.1e-3
    eps2 = 0.2e-5
    eps3 = 0.1e-6
    ncl = len( lelli )
    lthsv = (abs(sthetv - sthets) < eps2) and (abs(cphi - 1) < eps2)
    # lthsv: whether two directions coincide
    sumu  = 0
    sumd  = 0
    azs   = 0
    tgths = sthets/cthets
    tgthv = sthetv/cthetv
    for i_n in range(ncl): # end of loop at the end of the file (line 200)
        rlls2  = 0
        rllv2u = 0
        rllv2d = 0
        x8     = 0
        y8     = 0
        lellib = lelli[i_n]
        ulgi   = ulg[i_n]
        rcri   = rcr[i_n]
        htree  = htr[i_n]
        dbi    = dbh[i_n]
        hc1i   = hc1[i_n]
        hc2i   = hc2[i_n]
        shli   = shl[i_n]
        if lellib:
            cellb = hc1i*0.5
            hbase = htree - hc1i
        else:
            hbase = htree - (hc1i + hc2i)

        #  crown projection in the view direction (up)
        #  võra projektsioon vaatesuunas (üles)
        if z1 < htree:
            if lellib:
                szu1, vzui = pi11u (tgthv, z1, htree, cellb, rcri)
                # szu1: projected area
                #   vzui: volume of the upper part of the ellipsoid
            else:
                szu1, vzui = pi22u(tgthv, z1, htree, hc1i, hc2i, rcri)
            if szu1 > 0:
                rllv2u = vzui/szu1/cthetv
                azvu  = np.exp(-ulgi*rllv2u) # crown transmittance (down)
                                # or the gap probability inside the crown
                sl1   = ((htree + max(z1, hbase))*0.5 - z1)*tgthv # height of crown center relative to z1
            else:
                azvu  = 0
                sl1   = 0
        else:
            azvu  = 0
            sl1   = 0
            szu1  = 0
        #                        crown projection in the view direction (down)
        #                             võra projektsioon vaatesuunas (alla)
        if z3 > hbase:
            if lellib:
               szd3, vzdi = pi11d (tgthv, z3, htree, cellb, rcri)
            else:
               szd3, vzdi = pi22d(tgthv, z3, htree, hc1i, hc2i, rcri)
            if szd3 > 0:
               rllv2d = vzdi/szd3/cthetv
               azvd  = np.exp(-ulgi*rllv2d)
               sl3  = ((hbase + z3)*0.5 - hbase)*tgthv
            else:
               azvd  = 0
               sl3   = 0
        else:
            azvd  = 0
            sl3   = 0
            szd3  = 0
        #                      trunk projection in the view direction (up)
        #                          tu"ve projektsioon vaatesuunas (üles)
        slty1 = 0
        slty2 = 0
        if z1 < hbase:
            slty1 = (hbase - z1)*tgthv
            styx = stem(z1, hbase, dbi, htree)
            sty1  = styx*tgthv
            # sty1  = dbi*slty1
        else:
            sty1  = 0
        #                    trunk projection in the view direction (down)
        #                          tu"ve projektsioon vaatesuunas (alla)
        slty3 = min(z3, hbase)
        zx1 = 0
        zx2 = slty3
        styx = stem(zx1, zx2, dbi, htree)
        sty3  = styx*tgthv
        # sty3  = slty3*dbi
        #                              pooi üles binoomvalemiga, 28.05.1998
        #                                  esimene kordaja
        vaheg   = 1 - glmp[i_n]
        if abs(vaheg) > eps1:
            b11i = -np.log(1 - (1 -azvu)*vaheg)/vaheg
            # coefficient b1j in Eq. (7) of manual (ver. 05.2005)
        else:
            b11i = 1.0 - azvu

        if lthsv:
            # the two directions are the same, no need to calculate anything for the
            # Sun direction
            b2iz  = b11i
            b12i  = b11i
            szu2  = szu1
            scomm = szu1
            sty2  = sty1
            scty  = sty1
        else:
            #                   crown projection in the Sun direction
            #                         võra projektsioon Päikese suunas
            if z2 < htree:
                if lellib:
                    szu2, vzui = pi11u(tgths, z2, htree, cellb, rcri)
                else:
                  szu2, vzui = pi22u(tgths, z2, htree, hc1i, hc2i, rcri)
                if szu2 > 0:
                    rlls2 = vzui/szu2/cthets
                    azs   = np.exp(-ulgi*rlls2)
                    sl2   = ((htree + max(z2, hbase))*0.5 - z2)*tgths
                else:
                  azs   = 0
                  sl2   = 0
            else:
               azs  = 0
               sl2  = 0
               szu2 = 0
            #                    trunk projection in the Sun direction
            #                        tu"ve projektsioon thets suunas
            if z2 < hbase:
                slty2 = (hbase - z2)*tgths
                styx = stem(z2, hbase, dbi, htree)
                sty2  = styx*tgths
                # sty2  = slty2*dbi
            else:
               sty2  = 0

            #                           distance between projection centers
            #                        projektsioonitsentrite vaheline kaugus
            x8  = x2 + sl2
            y8  = y2
            x7  = x1 + sl1*cphi
            y7  = y1 + sl1*sphi
            sl8 = np.sqrt((x7 - x8)**2 + (y7 - y8)**2)
            #                            common part of trunk projections
            #                              tu"ve projektsioonide u"hisosa
            if (sty1 > 0) and (sty2 > 0):
                if cphi < -eps3:
                    scty = 0
                elif abs(y1 - y2) < eps1:
                    x5  = x2 + slty2
                    x6  = x1 + slty1*cphi
                    if abs(x1 - x2) < eps1:
                        tstphi = dbi/min(slty1, slty2)
                        if sphi > tstphi:
                            scty = dbi**2*0.5/sphi
                        elif (min(x5, x6) - max(x1, x2)) > 0:
                            scty = (min(x5, x6) - max(x1, x2))*dbi
                        else:
                            scty = 0
                    else:
                        if sphi > eps3:
                            scty = 0
                        elif (min(x5, x6) - max(x1, x2)) > 0:
                            scty = (min(x5, x6) - max(x1, x2))*dbi
                        else:
                            scty = 0
                        #
                    #
                else:
                  scty = 0
                #
            else:
                scty = 0
            #                            common part of crown projections
            #                      vo~ra projektsioonide u"hisosa pindala
            ss1 = max(szu1, szu2)
            ss2 = min(szu1, szu2)
            rc1 = np.sqrt(ss1/np.pi)
            rc2 = np.sqrt(ss2/np.pi)
            if sl8 > (rc1 + rc2):
                scomm = 0
            elif (sl8 < (rc1 - rc2)) or ( sl8 < eps2):
                scomm = ss2
            else:
                ppp   = (rc1 + rc2 + sl8)*0.5
                angl4 = np.arccos((rc1**2 + sl8**2 - rc2**2)/(2*rc1*sl8))
                angl8 = np.arccos((rc1**2 + rc2**2 - sl8**2)/(2*rc1*rc2))
                angl5 = (np.pi - angl4 - angl8)
                ss4   = angl4*rc1**2
                ss5   = angl5*rc2**2
                ss6   = np.sqrt(ppp*(ppp - rc1)*(ppp - rc2)*(ppp - sl8))
                scomm = ss4 + ss5 - 2*ss6
            #                        between-crown hot-spot correction
            #                           ! va"line hot-spoti korrektsioon
            rlls = rlls1 + rlls2
            rllv = rllv1 + rllv2u
            rl12 = (rllv - rlls)**2 + (1 - calph)*2*rlls*rllv
            if rl12 < 0:
                crr2 = 0
            else:
               crr2 = np.sqrt(rl12)/shli
            if crr2 < eps1:
                xtmp = 1 - crr2*0.5
            else:
                xtmp = rintegr(crr2)
            chs2 = np.exp((np.sqrt(rllv*rlls)*xtmp - chs1)*ulgi)

            azvs  = azvu*azs*chs2
            if abs(vaheg) > eps1:
               b12i = -np.log(1 - (1 - azs)*vaheg)/vaheg
               b2iz = -np.log(1 - (1 - azvs)*vaheg)/vaheg
            else:
               b12i   = 1 - azs
               b2iz   = 1 - azvs
        # end conditional (if lthsv)
        styx = min(sty1, sty2)
        styx = min(styx, scty)
        styx = sty1 + sty2 - styx
        sumu = sumu + stdns[i_n]*(b11i*(szu1 - scomm) + \
            b12i*(szu2 - scomm) + b2iz*scomm + styx)


        #                               pooi alla binoomvalemiga, 28.05.1998
        #                               esimene kordaja
        if abs(vaheg) > eps1:
            b11i = -np.log(1 - (1 -azvd)*vaheg)/vaheg
        else:
            b11i = 1 - azvd
        x7  = x3 - sl3*cphi
        y7  = y3 - sl3*sphi
        sl8 = np.sqrt((x7 - x8)**2 + (y7 - y8)**2)

        #                                tu"ve projektsioonide u"hisosa
        if (sty3 > 0) and (sty2 > 0):
            if -cphi < -eps3:
                scty = 0
            elif abs(y3 - y2) < eps1:
                x5  = x2 + slty2
                x6  = x3 - slty3*cphi
                if abs(x3 - x2) < eps1:
                    tstphi = dbi/min(slty3, slty2)
                    if sphi > tstphi:
                        scty = dbi**2*0.5/sphi
                    elif (min(x5, x6) - max(x3, x2)) > 0:
                        scty = (min(x5, x6) - max(x3, x2))*dbi
                    else:
                        scty = 0
                else:
                    if sphi > eps3:
                        scty = 0
                    elif (min(x5, x6) - max(x3, x2)) > 0:
                        scty = (min(x5, x6) - max(x3, x2))*dbi
                    else:
                        scty = 0
                #
            else:
                scty = 0
        else:
            scty = 0
        #                         vo~ra projektsioonide u"hisosa pindala
        ss1 = max(szd3, szu2)
        ss2 = min(szd3, szu2)
        rc1 = np.sqrt(ss1/np.pi)
        rc2 = np.sqrt(ss2/np.pi)
        if sl8 >= (rc1 + rc2):
            scomm = 0
        elif (sl8 <= (rc1 - rc2)) or (sl8 <= eps2):
            scomm = ss2
        else:
            ppp   = (rc1 + rc2 + sl8)*0.5
            s34   = ppp*(ppp - rc1)*(ppp - rc2)*(ppp - sl8)
            sl4   = 2/sl8*np.sqrt(s34)
            angl1 = 2*np.arcsin(sl4/rc1)
            angl2 = 2*np.arcsin( min(sl4/rc2,1) ) # the argument can exceed 1 and produce nan
            ss4   = rc1**2*0.5*(angl1 - np.sin(angl1))
            ss5   = rc2**2*0.5*(angl2 - np.sin(angl2))
            scomm = ss4 + ss5
        #                                 va"line hot-spoti korrektsioon
        rlls = rlls1 + rlls2
        rllv = rllv3 + rllv2d
        rl12 = (rllv - rlls)**2 + (1 + calph)*2*rlls*rllv
        if rl12 <= 0:
            crr2 = 0
        else:
            crr2 = np.sqrt(rl12)/shli
        if crr2 < eps1:
            xtmp = 1.0 - crr2*0.5
        else:
            xtmp = rintegr(crr2)
        chs2 = np.exp((np.sqrt(rllv*rlls)*xtmp - chs3)*ulgi)

        azvs  = azvd*azs*chs2
        if abs(vaheg) > eps1:
            b12i = -np.log(1 - (1- azs)*vaheg)/vaheg
            b2iz = -np.log(1 - (1 - azvs)*vaheg)/vaheg
        else:
            b12i   = 1 - azs
            b2iz   = 1 - azvs

        styx = min(sty3, sty2)
        styx = min(styx, scty)
        styx = sty3 + sty2 - styx
        sumd = sumd + stdns[i_n]*(b11i*(szd3 - scomm) + \
            b12i*(szu2 - scomm) + b2iz*scomm + styx)

    pooui = np.exp(-sumu)
    poodi = np.exp(-sumd)
    return pooui, poodi
#   -------------------------- end spooi


def bck3(lelli, icl, ul, shl, stdns, htr, dbh, hc1, hc2, rcr,
    ulg, glmp, sthets, cthets, sthetv, cthetv, sphi,
    cphi, calph, xi, yj, zk):
    """Bidirectional gap probability in a tree crown for the spherical orientation of leaves.

    subroutine bck3 Jo~eveer, (mai 88) Andrese progr.-st ilips xi, yj, zk
    A. Kuusk  4.01.1986. Muudetud 11.05.1998. Pyhton translation Matti Mõttus 2022
    1. ja"rku hajumise heleduskoefitsient (vo~rad)
    A. Kuusk, The hot spot effect in plant canopy reflectance.
    In: R. Myneni and J. Ross (Eds.), Photon-Vegetation Interactions.
    Applications in Optical Remote sensing and Plant Ecology.
    Springer, Berlin, 1991, 139-159.

    Args:
    lelli: Flag. True -> ellipsoid; False -> -- koonus + silinder
    calph: cos(alph), alph - nurk kahe kiire vahel
    shl: hot spot parameter, lehe (kasvu) suurus (m)
    sthet: sin(thet?)
    cthet: cos(thet?)
    sphi: sin(phi)
    cphi: cos(phi), phi on nurk kahte kiirt la"bivate vertikaaltasandite vahel.
    rcr: aell = võra raadius
    hcl: vhelko = võra ko~rgus
    htr: puu ko~rgus

    Returns:
    sxui: gap probability, view from above
    sxdi: gap probability, view from below
    """
    vhelko = hc1[icl]
    lellib = lelli[icl]
    aell   = rcr[icl]
    htree  = htr[icl]
    shli   = shl[icl]
    if lellib:
        cell   = 0.5*vhelko
        rllvu = rlips(xi, yj, zk, sthetv, cthetv, sphi, cphi, aell, cell)
        rllvd = rlips(xi, yj, zk, sthetv, -cthetv, sphi, -cphi, aell, cell)
        rlls = rlips (xi, yj, zk, sthets, cthets, 0.0, 1.0, aell, cell)
    else:
        cell   = vhelko
        rllvu = rcone(xi, yj, zk, sthetv, cthetv, sphi, cphi, aell, vhelko, hc2[icl])
        rllvd = rcone(xi, yj, zk, sthetv, -cthetv, sphi, -cphi,aell, vhelko, hc2[icl])
        rlls = rcone(xi, yj, zk, sthets, cthets, 0.0, 1, aell, vhelko, hc2[icl])
    #   hot spot correction, hot-spoti korrektsioon
    rl12u = (rllvu - rlls)**2 + (1 - calph)*2*rlls*rllvu
    if rl12u < 0:
        crr0 = 0
    else:
        crr0 = np.sqrt(rl12u)/shli

    if crr0 < 0.001:
        xtmp = 1 - crr0*0.5
    else:
        xtmp = rintegr(crr0)
    chs1 = np.sqrt(rllvu*rlls)*xtmp
    xx  = (chs1 - rlls - rllvu)*ul*0.5
    bgpu = np.exp(xx)

    rl12d = (rllvd - rlls)**2 + (1 + calph)*2*rlls*rllvd
    if rl12d < 0:
        crr0 = 0
    else:
        crr0 = np.sqrt(rl12d)/shli

    if crr0 < 0.001:
        xtmp = 1.0 - crr0*0.5
    else:
        xtmp = rintegr(crr0)

    chs3 = np.sqrt(rllvd*rlls)*xtmp
    xx  = (chs3 - rlls - rllvd)*ul*0.5
    bgpd = np.exp(xx)
    #                bgp - kahe suuna ava to~ena"osus võra sees
    #                      bidirectional gap probability
    #  To~ena"osuste pooi arvutamine.  (A. Jo~eveer)
    #  Vaba vaatesuuna tõenäosused väljaspool võra korraga päikese suunas
    #  ja vaatesuunas nii üles kui alla, pooui ja poodi, vastavalt.

    x1  = xi + rllvu*sthetv*cphi
    y1  = yj + rllvu*sthetv*sphi
    z1  = zk + rllvu*cthetv + htree - cell
    #                htree - puu ko~rgus,
    #                (htree - cell) - võra koordinaadistiku alguspunkt
    x2  = xi + rlls*sthets
    y2  = yj
    z2  = zk + rlls*cthets + htree - cell

    x3  = xi - rllvd*sthetv*cphi
    y3  = yj + rllvd*sthetv*sphi
    z3  = zk - rllvd*cthetv + htree - cell

    pooui, poodi = spooi(lelli, ulg, shl, stdns, htr, hc1, hc2, rcr, dbh,
        glmp, x1, y1, z1, x2, y2, z2, x3, y3, z3, rlls, rllvu, rllvd,
        sthets, cthets, sthetv, cthetv, sphi, cphi, calph, chs1, chs3)
    sxui  = bgpu*pooui
    sxdi  = bgpd*poodi
    return sxui, sxdi
#   -------------------------- end bck3


def stem(z1, z2, dbh, htree):
    """ stem vertical cross-section area

    subroutine stem  A. Kuusk 20.04.2005, python translation Matti Mõttus 2022
    Tüve pikilõike pindala. Tüve moodustaja valemid on pärit lätlaste 1988. a.
    normatiivdokumentidest. Pindala arvutatakse nivoost z1 nivooni z2 ristlõike
    diameetri integreerimisega.

    Args:
    dbh: rinnasdiameeter
    htree: puu kõrgus

    Returns:
    styvi
    """
    #  need parameetrid on männi jaoks
    #  parameters for Scots pine
    ati = [ 118.981, -277.578, 1140.525, -3037.487, 4419.682, -3361.78, 997.657 ]
    ht0 = 26.0,
    dt0 = 30.0
    pty = 0.007
    qty = -0.007

    xz0 = 1.3/htree
    xz1 = z1/htree
    xz2 = z2/htree
    eet = pty*(htree - ht0) + qty*(dbh - dt0)

    sum1 = 0.0
    for j in range(7):
        sum1 += ati[j]*(xz0**j)
    f13 = sum1*(1.0 + eet*(xz0*xz0 - 0.01))

    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    sum4 = 0.0
    for  j in range(1,8):
        sum1 = sum1 + ati[j-1]/j*(xz1**j)
        sum2 = sum2 + ati[j-1]/j*(xz2**j)
        sum3 = sum3 + ati[j-1]/(j+3)*(xz1**(j + 3))
        sum4 = sum4 + ati[j-1]/(j+3)*(xz2**(j + 3))

    styvi = (sum2 - sum1)*(1.0 - 0.01*eet)
    styvi = styvi + eet*(sum4 - sum3)
    styvi = styvi*htree*dbh/f13
    return styvi
#   -------------------------- end stem


def rlips (xi, yj, zk, sthrl, cthrl, sphrl, cphrl, aell, cell):
    """ distance from (xi,yi,zk) to the surface of the ellipsoid in the direction thet, phi

    subroutine rlips, A. Kuusk   29.12.1985; python translation Matti Mõttus (2022)
    kaugus p - tist xi, yj, zk ellipsoidi pinnani suunas thet, phi;
    phi on nurk x-telje ja kiirt (thet, phi) la"biva vertikaaltasandi vahel.
    Koordinaatide algus on ellipsoidi tsentris.
    Origin of coordinates is in the center of the ellipsoid of rotation

    Args:
    sthrl = sin(thet)
    cthrl = cos(thet), 0 <= thet <= pi/2
    sphrl = sin(phi),
    cphrl = cos(phi),  0 <= phi  <= pi
    aell: horizontal semi-axis of the ellipoid of rotation (crown radius)
    cell: vertical semi-axis of the ellipoid of rotation (1/2 crown length9

    Returns:
    rlout
    """
    eps = 0.1e-4

    a2    = aell*aell
    c2    = cell*cell
    aaa   = c2*(sthrl**2) + a2*(cthrl**2)
    bbb   = 2.0*(c2*sthrl*(xi*cphrl + yj*sphrl) + a2*cthrl*zk)
    ccc   = c2*(xi**2 + yj**2) + a2*(zk**2 - c2)
    det   = (bbb**2 - 4.0*aaa*ccc)
    if abs(det) < eps:
        rlout = -bbb*0.5/aaa
    else:
        rlout = (np.sqrt(bbb**2 - 4.0*aaa*ccc) - bbb)*0.5/aaa

    if rlout < 0:
            rlout = 0
    return rlout
#   --------------------------- end rlips


def rcone(xi, yj, zk, sthrl, cthrl, sphrl, cphrl, aell, vhelko, vhcyl):
    """ Distance from the point (xi, yj, zk) to the surface of the crown

    subroutine rcone A. Kuusk, 7.05.1998
    kaugus p-tist (xi, yj, zk) võra pinnani suunas thet, phi

    Args:
    vhelko - koonilise osa ko~rgus,  aell - võra raadius,
    phi on nurk x-teljest, sphrl = sin(phi), cphrl = cos(phi)
    sthrl = sin(thet), cthrl = cos(thet)

    Returns:
    rlout
    """
    eps = 0.1e-6
    eps2 = 0.001

    z3 = vhelko + eps2
    z0 = zk - eps2
    if (cthrl >= 0):
        # upper hemisphere, ülemine poolsfäär
        if (sthrl < eps):
            # ray is vertical, kiir on vertikaalne
            rlout = vhelko - vhelko/aell*np.sqrt(xi**2 + yj**2) - zk
        else:
            ae2 = aell**2
            hh2 = vhelko**2
            if (zk >= 0):
                # (x,y,z) is in the cone / on koonuses
                aaa = ae2/hh2*(cthrl**2) - sthrl**2
                bbb = 2.0*ae2*cthrl/vhelko*(zk/vhelko - 1.0) - 2.0*sthrl*(xi*cphrl + yj*sphrl)
                ccc = (aell*(1.0 - zk/vhelko))**2 - (xi**2 + yj**2)
                det = bbb**2 - 4.0*aaa*ccc
                if (abs(det) < eps):
                    rlout = -bbb*0.5/aaa
                else:
                    rlout = (np.sqrt(det) - bbb)*0.5/aaa
                    z2    = rlout*cthrl + zk
                    if (z2 < z0 or (z2 > z3) ):
                        rlout = (-np.sqrt(det) - bbb)*0.5/aaa
            else:
                # (x,y,z) is in the cylinder, on silindris
                rl1 = 0.0
                aaa = sthrl**2
                bbb = 2.0*sthrl*(xi*cphrl + yj*sphrl)
                ccc = xi**2 + yj**2 - ae2
                det = bbb**2 - 4.0*aaa*ccc
                if (abs(det) < eps):
                    rl1 = -bbb*0.50/aaa
                else:
                    rl1 = (np.sqrt(det) - bbb)*0.50/aaa
                    z2  = rl1*cthrl + zk
                    if (z2 < z0):
                        rl1 = (-np.sqrt(det) - bbb)*0.5/aaa
                z2  = rl1*cthrl + zk
                if (z2 <= 0.0):
                    # ray exits through side of cylinder, kiir väljub läbi silindri külgpinna
                    rlout = rl1
                else:
                    # ray exits through cone, kiir väljub läbi koonuse
                    rl2 = 0.0
                    rl1 = abs(zk)/cthrl
                    x3  = xi + rl1*sthrl*cphrl
                    y3  = yj + rl1*sthrl*sphrl
                    aaa = ae2/hh2*(cthrl**2) - sthrl**2
                    bbb = -2.0*ae2*cthrl/vhelko - 2.0*sthrl*(x3*cphrl + y3*sphrl)
                    ccc = ae2 - (x3**2 + y3**2)
                    det = bbb**2 - 4.0*aaa*ccc
                    if (abs(det) < eps):
                        rl2 = -bbb*0.5/aaa
                    else:
                        rl2 = (np.sqrt(det) - bbb)*0.5/aaa
                        z2  = rl2*cthrl
                        if (z2 < 0.0 or (z2 > z3) ):
                            rl2 = (-np.sqrt(det) - bbb)*0.5/aaa
                    rlout = rl1 + rl2
                #
            #
        #
    else:
        # lower hemisphere, alumine poolsfäär
        if (sthrl < eps):
            # kiir on vertikaalne
            rlout = vhcyl + zk
        else:
            ae2 = aell**2 + eps
            if (zk >= 0.0):
                # (x,y,z) is in the cone, on koonuses
                rl1 = -zk/cthrl
                x2  = xi + rl1*sthrl*cphrl
                y2  = yj + rl1*sthrl*sphrl
                if ((x2**2 + y2**2) >= ae2):
                    # ray exits through the surface of the cone, kiir väljub läbi koonuse pinna
                    hh2 = vhelko**2
                    aaa = aell**2/hh2*(cthrl**2) - sthrl**2
                    bbb = 2.0*aell**2*cthrl/vhelko*(zk/vhelko - 1.0) - 2.0*sthrl*(xi*cphrl + yj*sphrl)
                    ccc = (aell*(1.0 - zk/vhelko))**2 - (xi**2 + yj**2)
                    det = bbb**2 - 4.0*aaa*ccc
                    if (abs(det) < eps):
                        rlout = -bbb*0.5/aaa
                    else:
                        rlout = (np.sqrt(det) - bbb)*0.5/aaa
                        if (rlout < 0.0):
                            rlout = (-np.sqrt(det) - bbb)*0.5/aaa

                else:
                    # ray enters the cylinder, kiir siseneb silindrisse
                    rl1 = -(vhcyl + zk)/cthrl
                    x2  = xi + rl1*sthrl*cphrl
                    y2  = yj + rl1*sthrl*sphrl
                    if ((x2**2 + y2**2) <= ae2):
                        # ray exits through cylinder base, kiir väljub läbi silindri põhja
                        rlout = rl1
                    else:
                        # ray exits through cylinder wall, kiir väljub läbi silindri seina
                        aaa = sthrl**2
                        bbb = 2.0*sthrl*(xi*cphrl + yj*sphrl)
                        ccc = xi**2 + yj**2 - aell**2
                        det = bbb**2 - 4.0*aaa*ccc
                        if (abs(det) < eps):
                            rl1 = -bbb*0.5/aaa
                        else:
                            rl1 = (np.sqrt(det) - bbb)*0.5/aaa
                        z2  = rl1*cthrl + zk
                        if (z2 > 0.0):
                            rl1 = (-np.sqrt(det) - bbb)*0.5/aaa
                        rlout = rl1
                    #
                #
            else:
                # (x,y,z) is in the cylinder, on silindris
                rl1 = -(vhcyl + zk)/cthrl
                x2  = xi + rl1*sthrl*cphrl
                y2  = yj + rl1*sthrl*sphrl
                if ((x2**2 + y2**2) <= ae2):
                    # ray exits trhough sylinder bottom, kiir väljub läbi silindri põhja
                    rlout = rl1
                else:
                    # ray exits through sylinder wall, kiir väljub läbi silindri seina
                    aaa = sthrl**2
                    bbb = 2.0*sthrl*(xi*cphrl + yj*sphrl)
                    ccc = xi**2 + yj**2 - aell**2
                    det = bbb**2 - 4.0*aaa*ccc
                    if (abs(det) < eps):
                        rl1 = -bbb*0.5/aaa
                    else:
                        rl1 = (np.sqrt(det) - bbb)*0.5/aaa
                    z2  = rl1*cthrl + zk
                    if (z2 > zk + eps):
                        rl1 = (-np.sqrt(det) - bbb)*0.5/aaa
                    rlout = rl1
                #
            #
        #
    #
    if (rlout < 0.0):
        rlout = 0.0
    return rlout
#   --------------------------- end rcone
