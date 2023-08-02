# helper functions used by the frt class
#  these functions are not dependent on any specifics of frt and maybe useful elsewhere, too
import numpy as np
import os
import sys

def rintegr(xxx):
    return (1 - np.exp(-xxx))/xxx

def eint( xy ):
    return np.where( xy==0, 1, (1.0 - np.exp(-xy))/xy)
# eint will create an invalid value warning if xy==0 (view and solar angles coincide)

def CI_minfun( FGI, clmpst ):
    """ minimization function for finding Fischer's grouping index FGI from
        the crown-level clumping index clmpst
    """
    if FGI == 1:
        return clmpst-1
    else:
        return clmpst+np.log(FGI)/(1-FGI)
#   -------------------------- end CI_minfun


def gauleg_py(x1, x2, n):
    """Gauss-Legendre n-point quadrature. Press et al., 1992, Chpt. 4.5

    WARNING: use gauleg_numpy() instead!

    converted from subroutine gauleg(x1, x2, x, w, n)
    original frt.f: call gauleg(x1, x2, zctst, acztst, ntstz)

    Args:
    in: x1, x2 interval start and end
    n: number of nodes -- must be even!

    Returns:
    x: points
    w: weights (guess)
    """

    x = np.empty( n )
    w = np.empty_like( x )
    eps = 1.0e-10
    m  = int( np.ceil((n + 1)/2) )
    xm = 0.5*(x2 + x1)
    xl = 0.5*(x2 - x1)
    for i in range(m):
        z = np.cos(np.pi*( i + 0.75 )/( n + 0.5 ))
        while True: # until loop in fortran
            p1 = 1
            p2 = 0
            for j in range(n):
                p3 = p2
                p2 = p1
                p1 = ( (2*j + 1)*z*p2 - j*p3 )/(j+1)
            pp = n*(z*p1 - p2)/(z*z - 1)
            z1 = z
            z  = z1 - p1/pp
            if (abs(z - z1) < eps):
                break
        if (abs(z) < eps):
            z = 0
        x[i]     = xm - xl*z
        x[n-i-1] = xm + xl*z
        w[i]     = 2*xl/( (1 - z*z)*pp*pp )
        w[n-i-1] = w[i]
    return x, w
#   -------------------------- end gauleg_py


def gauleg_numpy(x1,x2,n):
    """Gauss-Legendre quadrature with n points between x1 and x2
    Same as gauleg_py() above but with numpy
    """
    x,w = np.polynomial.legendre.leggauss(n)
    x = x1 + (x+1)*(x2-x1)/2
    w *= (x2-x1)/2
    return x, w
#   -------------------------- end gauleg_numpy

def cubell9_py( ):
    """ knots and weights of the cubature in a sphere

    Converted from fortran77 frt code written by A. Kuusk  08.01.2003
    Mysovskih,I.P. Interpolyatsionnye kubaturnye formuly,
    Nauka, Moskva 1981, Chpt. 5, Algorithm 16.31.

    Args: no inputs?

    Returns:
    ntst: cubature size
    xtst, ytst, ztst: cubature knots
    atst: weights?

    """

    i1 = [ 1,  1,  1,  1,  1,  2,  2,  2,  2,  2,
           3,  3,  3,  4,  4,  5,  5,  6,  6,  7 ]
    i2 = [ 3,  3,  4,  5,  6,  8,  8,  9, 10, 11,
           4,  7,  8,  5,  8,  6,  9,  7, 10, 11 ]
    i3 = [ 4,  7,  5,  6,  7,  9, 12, 10, 11, 12,
           8, 12, 12,  9,  9, 10, 10, 11, 11, 12 ]

    nface = 20
    ntst = 45
    vol  = np.pi*4/3
    atst = [None]*ntst
    atst[0] = vol*2096/42525
    atst[1:13] = [ vol*(491691 + 54101*np.sqrt(31))/21.0924e6 ]*12
    atst[13:25] = [ vol*(491691 - 54101*np.sqrt(31.))/21.0924e6 ]*12
    atst[25:45] = [ vol*1331/68.04e3 ]*20

    alph = np.sqrt(81 - 6*np.sqrt(31))/11
    beta = np.sqrt(81 + 6*np.sqrt(31))/11
    gamm = 3/np.sqrt(11)
    #                                   *****    center of sphere
    xtst = [None]*ntst
    xtst[0] = 0
    ytst = [None]*ntst
    ytst[0] = 0
    ztst = [None]*ntst
    ztst[0] = 0
    #                                   *****    verteces of icosahedron (12)
    xi = [None]*ntst
    xi[1] = 0
    xi[2] = 0
    yi = [None]*ntst
    yi[1] = 0
    yi[2] = 0
    zi = [None]*ntst
    zi[1] = 1
    zi[2] = -1

    for i in range(3,8):
        xi[i] = np.cos( (i-3)*2*np.pi /5 )*2/np.sqrt(5)
        xi[i+5] = np.cos( (2*(i-3)+1 )*2*np.pi/5)*2/np.sqrt(5)
        yi[i] = np.sin( (i-3)*2*np.pi/5 )*2/np.sqrt(5)
        yi[i+5] = np.sin( (2*(i-3)+ 1)*2*np.pi/5)*2/np.sqrt(5)
        zi[i] = 1/np.sqrt(5)
        zi[i+5] = -1/np.sqrt(5)

    for i in range(1,13):
        xtst[i] = xi[i]*alph
        ytst[i] = yi[i]*alph
        ztst[i] = zi[i]*alph
        xtst[i+12] = xi[i]*beta
        ytst[i+12] = yi[i]*beta
        ztst[i+12] = zi[i]*beta
    #                                    ***** projections of facets' centers
    i  = 0
    xz = xi[ i1[i] ] + xi[ i2[i] ] + xi[ i3[i] ]
    yz = yi[ i1[i] ] + yi[ i2[i] ] + yi[ i3[i] ]
    zz = zi[ i1[i] ] + zi[ i2[i] ] + zi[ i3[i] ]
    rz = np.sqrt(xz**2 + yz**2 + zz**2)
    xtst[25] = xz/rz*gamm
    ytst[25] = yz/rz*gamm
    ztst[25] = zz/rz*gamm
    for i in range (1,nface):
        xtst[i+25] = (xi[ i1[i] ] + xi[ i2[i] ] + xi[ i3[i] ])/rz*gamm
        ytst[i+25] = (yi[ i1[i] ] + yi[ i2[i] ] + yi[ i3[i] ])/rz*gamm
        ztst[i+25] = (zi[ i1[i] ] + zi[ i2[i] ] + zi[ i3[i] ])/rz*gamm
    return ntst, xtst, ytst, ztst, atst
#   -------------------------- end cubell9_py


def cubell9_py_fixed_DONOTUSE( ):
    """ knots and weights of the cubature in a sphere
    XXX -- needs to be re-run, bug fixed-- XXX
    The numerical output from cubell9_py(). If the output is constant, why not just tabulate it?
    """
    ntst = 45
    xtst = [0.0,                  0.0,                  0.0,
            0.5609520468185168,   0.17334371549633287, -0.4538197389055912,
        -0.45381973890559124,  0.17334371549633273,  0.17334371549633287,
        -0.45381973890559124,  0.5609520468185168,  -0.4538197389055911,
            0.1733437154963326,   0.0,                  0.0,
            0.8697167247646821,   0.2687572482444055,  -0.7036156106267464,
        -0.7036156106267465,   0.2687572482444053,   0.2687572482444055,
        -0.7036156106267465,   0.8697167247646821,  -0.7036156106267463,
            0.26875724824440506,  0.444237896264228,    0.444237896264228,
        -0.16968377728218514, -0.5491082379640857,  -0.16968377728218526,
        -0.16968377728218514,  0.20974068339971527,  0.0648134355823275,
            0.06481343558232758, -0.16968377728218526,  0.5491082379640858,
            0.5491082379640855,   0.5491082379640857,  -0.444237896264228,
        -0.06481343558232745, -0.20974068339971538, -0.20974068339971538,
        -0.44423789626422805, -0.20974068339971533, -0.06481343558232762]

    ytst = [0.0,                  0.0,                  0.0,
            0.0,                  0.5334970994558544,   0.32971934036320116,
        -0.32971934036320105, -0.5334970994558546,   0.5334970994558544,
        -0.32971934036320105,  0.0,                  0.3297193403632012,
        -0.5334970994558546,   0.0,                  0.0,
            0.0,                  0.8271497584183294,   0.5112066644887922,
        -0.5112066644887919,  -0.8271497584183295,   0.8271497584183294,
        -0.5112066644887919,   0.0,                  0.5112066644887923,
        -0.8271497584183296,  -0.3227577241875956,  -0.3227577241875956,
            0.5222329678670937,   0.0,                 -0.5222329678670936,
            0.1232824805080975,   0.0,                 -0.19947524367949812,
            0.1994752436794981,  -0.12328248050809744,  0.645515448375191,
        -0.6455154483751913,   0.0,                  0.32275772418759563,
            0.446040204695693,    0.0,                  0.0,
        -0.3227577241875954,   0.0,                 -0.446040204695693]

    ztst = [0.0,                  0.6271634544019241,  -0.6271634544019241,
            0.2804760234092584,   0.2804760234092584,   0.2804760234092584,
            0.2804760234092584,   0.2804760234092584,  -0.2804760234092584,
        -0.2804760234092584,  -0.2804760234092584,  -0.2804760234092584,
        -0.2804760234092584,   0.972372858871152,   -0.972372858871152,
            0.43485836238234105,  0.43485836238234105,  0.43485836238234105,
            0.43485836238234105,  0.43485836238234105, -0.43485836238234105,
        -0.43485836238234105, -0.43485836238234105, -0.43485836238234105,
        -0.43485836238234105,  0.7187920152462709,   0.7187920152462709,
            0.7187920152462709,   0.7187920152462709,   0.7187920152462709,
        -0.7187920152462709,  -0.7187920152462709,  -0.7187920152462709,
        -0.7187920152462709,  -0.7187920152462709,   0.16968377728218517,
            0.16968377728218517, -0.16968377728218517,  0.16968377728218517,
        -0.16968377728218517,  0.16968377728218517, -0.16968377728218517,
            0.16968377728218517, -0.16968377728218517, -0.16968377728218517]

    atst =  [0.20645982996430978, 0.1574664151562833,   0.1574664151562833,
            0.1574664151562833,   0.1574664151562833,   0.1574664151562833,
            0.1574664151562833,   0.1574664151562833,   0.1574664151562833,
            0.1574664151562833,   0.1574664151562833,   0.1574664151562833,
            0.1574664151562833,   0.03782577014094462,  0.03782577014094462,
            0.03782577014094462,  0.03782577014094462,  0.03782577014094462,
            0.03782577014094462,  0.03782577014094462,  0.03782577014094462,
            0.03782577014094462,  0.03782577014094462,  0.03782577014094462,
            0.03782577014094462,  0.0819412075627673,   0.0819412075627673,
            0.0819412075627673,   0.0819412075627673,   0.0819412075627673,
            0.0819412075627673,   0.0819412075627673,   0.0819412075627673,
            0.0819412075627673,   0.0819412075627673,   0.0819412075627673,
            0.0819412075627673,   0.0819412075627673,   0.0819412075627673,
            0.0819412075627673,   0.0819412075627673,   0.0819412075627673,
            0.0819412075627673,   0.0819412075627673,   0.0819412075627673]
    return ntst, xtst, ytst, ztst, atst
#   -------------------------- end cubell9_py_fixed


def cubcirc_py():
    """knots and weights of the quadrature in a circle

    Converted from fortran77 frt code written by A. Kuusk  27.08. 2002
    Vysotskikh, I.P., Nauka, M., 1981, Ch. 5, P-t. 16.33.

    Returns:
    ntst: quadrature size
    xtst, ytst: quadrature knots
    atst: quadrature weights
    """

    ri = [ 1.0, 2.0, 4.0, 5.0, 7.0, 8.0, 10.0, 11.0 ]

    ntst = 28

    r1 = np.sqrt((10.0 - np.sqrt(10.0))/15.0)
    r2 = np.sqrt((10.0 + np.sqrt(10.0))/15.0)
    r3 = np.sqrt((31.0 - np.sqrt(601.0))/60.0)
    r4 = np.sqrt(3.0/5.0)
    r5 = np.sqrt((31.0 + np.sqrt(601.0))/60.0)
    c1 = (340.0 + 25.0*np.sqrt(10.0))*np.pi/10368.0
    c2 = (340.0 - 25.0*np.sqrt(10.0))*np.pi/10368.0
    c3 = (857.0 + 12707.0/np.sqrt(601.0))*np.pi/20736.0
    c4 = 125.0*np.pi/3456.0
    c5 = (857.0 - 12707.0/np.sqrt(601.0))*np.pi/20736.0

    xtst = [None]*ntst
    ytst = [None]*ntst
    atst = [None]*ntst
    for k in range(8):
        xtst[k]   = r1*np.cos(np.pi*ri[k]/6.0)
        ytst[k]   = r1*np.sin(np.pi*ri[k]/6.0)
        atst[k]   = c1
        xtst[8+k] = r2*np.cos(np.pi*ri[k]/6.0)
        ytst[8+k] = r2*np.sin(np.pi*ri[k]/6.0)
        atst[8+k] = c2

    for k in range(4):
         xtst[16+k] = r3*np.cos(np.pi*(k+1)/2.0)
         ytst[16+k] = r3*np.sin(np.pi*(k+1)/2.0)
         atst[16+k] = c3
         xtst[20+k] = r4*np.cos(np.pi*(k+1)/2.0)
         ytst[20+k] = r4*np.sin(np.pi*(k+1)/2.0)
         atst[20+k] = c4
         xtst[24+k] = r5*np.cos(np.pi*(k+1)/2.0)
         ytst[24+k] = r5*np.sin(np.pi*(k+1)/2.0)
         atst[24+k] = c5

    return ntst, xtst, ytst, atst
#   --------------------------- end cubcirc_py


def cubcirc_py_fixed():
    """knots and weights of the quadrature in a circle

    The numerical output from cubcirc_py(). If the output is constant, why not just tabulate it?
    """
    ntst = 28

    xtst= [  0.584710284663765,   0.33758264024856743, -0.3375826402485672,
            -0.584710284663765,  -0.5847102846637651,  -0.33758264024856766,
             0.33758264024856743, 0.5847102846637647,   0.811242185175561,
             0.46837089398909043,-0.4683708939890901,  -0.811242185175561,
            -0.811242185175561,  -0.4683708939890907,   0.46837089398909043,
             0.8112421851755607,  0.0,                 -0.32875265919678565,
             0.0,                 0.32875265919678565,  0.0,
            -0.7745966692414834,  0.0,                  0.7745966692414834,
             0.0,                -0.9619017737816972,   0.0,
             0.9619017737816972 ]

    ytst = [ 0.3375826402485673,  0.5847102846637648,   0.584710284663765,
             0.3375826402485673, -0.3375826402485672,  -0.5847102846637648,
            -0.5847102846637648, -0.33758264024856766,  0.46837089398909026,
             0.8112421851755609,  0.811242185175561,    0.46837089398909026,
            -0.46837089398909004,-0.8112421851755608,  -0.8112421851755609,
            -0.4683708939890907,  0.32875265919678565,  0.0,
            -0.32875265919678565, 0.0,                  0.7745966692414834,
             0.0,                -0.7745966692414834,   0.0,
             0.9619017737816972,  0.0,                 -0.9619017737816972,
             0.0 ]

    atst = [ 0.12697783650322456, 0.12697783650322456,  0.12697783650322456,
             0.12697783650322456, 0.12697783650322456,  0.12697783650322456,
             0.12697783650322456, 0.12697783650322456,  0.07906797796832823,
             0.07906797796832823, 0.07906797796832823,  0.07906797796832823,
             0.07906797796832823, 0.07906797796832823,  0.07906797796832823,
             0.07906797796832823, 0.2083682752319387,   0.2083682752319387,
             0.2083682752319387,  0.2083682752319387,   0.11362820651004749,
             0.11362820651004749, 0.11362820651004749,  0.11362820651004749,
             0.0513100527123565,  0.0513100527123565,   0.0513100527123565,
             0.0513100527123565 ]
    return ntst, xtst, ytst, atst
#   --------------------------- end cubcurc_py_fixed


def hetk8o( thets, thetv_vec, phi_vec, stdns, uuu, bdgfu, bdgfd, btr1uk, btr1dk,
    pkhair, RefrIndex, TrunkRefl, LeafBranchReflLamb, LeafBranchTrans ):
    """Single scattering of crowns as viewed from above and below, f77 hetk8o

    võrade ühekordse hajumise heledus alt ja ülalt vaadates A. Kuusk  27.11.2000
    [ Eq. (3) of manual ver. 09.2002 ]. Python translation Matti Mõttus 2022
    Compute the scattering coefficients and apply them to the integrated gap fractions
    thetv_vec and phi_vec are vectors, all directions ( thetv_vec x phi_vec ) are sampled

    Args:
    thets: solar zenith angle [rad]
    thetv_vec: vector of view zenith angles [rad], e.g. knots of G-L quadrature (zenith angle)
    phi_vec: vector of relative view azimuth angles for the quadrature [rad]    thets: sun zenith angle [rad]
    stdns: stand density for each tree class
    uuu: leaf+branch area density for each tree class
    bdgfu, bdgfd, btr1uk, btr1dk: gap fractions integrated over the crowns of each class
    pkhair: correction to fresnel reflectance (?)
    RefrIndex: wax refractive index (for each tree class, f77 rnlf)
    TrunkRefl: trunk reflectance factor (for each tree class, f77 rtrnk)
    LeafBranchReflLamb: leaf (+brance) Lambertian reflectance factor (for each tree class, f77 rrs)
    LeafBranchTrans: leaf transmittance factor (for each tree class, f77 ttt)

    Returns:
    bi0u: single scattering reflectance component (tree crowns), 3D ndarray, thetv_vec x phi_vec X spectrum
    bi0d: single scattering transmittance component (tree crowns), 3D ndarray, thetv_vec x phi_vec X spectrum
    """

    eps = 0.01
    sthetv  = np.sin(thetv_vec)
    cthetv  = np.cos(thetv_vec)
    cthets  = np.cos(thets)
    sthets  = np.sin(thets)
    tgths   = sthets/cthets

    ncl = len(stdns)
    nth = len(thetv_vec)
    nphi = len(phi_vec)
    nwl = len(LeafBranchReflLamb[0])
    # initilize output variables
    bi0u = np.zeros( (nth,nphi,nwl) )
    bi0d = np.zeros( (nth,nphi,nwl) )
    # Vo~rade u"hekordse hajumise heleduse summeerimine u"le suurusklasside
    # Single-scattering reflectance [ Eq. (3) of manual ver. 09.2002]
    for iphi,phi in enumerate(phi_vec):
        cphi    = np.cos(phi)
        sphi    = np.sin(phi)
        calph   = sthetv*sthets*cphi + cthetv*cthets
        alpha   = np.arccos(calph)
        alpha[ np.nonzero( np.abs(alpha)<eps ) ] = 0
        alpha[ np.nonzero( np.abs(alpha)<0 ) ] += np.pi
        salph   = np.sin(alpha)
        calp2a  = np.cos(0.5*alpha) # cosine of incidence angle for leaves which cause Fresnel erflection into the sensor
        calp2b  = np.sqrt( (1 - calph)/2 )
        gammr = (salph + (np.pi - alpha)*calph)/(3*np.pi)
        gammt = (salph - alpha*calph)/(3*np.pi)
        for icl in range(ncl):
            # ------ upper hemisphere
            # calculate single-scattering reflectance from leaves
            # Gamma,  spherical leaf orientation, bi-Lambertian leaves, upper hemisphere
            gammd = gammr[..., np.newaxis]*LeafBranchReflLamb[icl] + gammt[...,np.newaxis]*LeafBranchTrans[icl]
            # gammd is 2D ndarray ( 1 x nwl )
            # add Fresnel reflectance
            rfr = gmfres( calp2a[..., np.newaxis], RefrIndex[icl], pkhair )

            gammout = gammd + rfr*0.125 # 1/8  comes from e.g. Kuusk, 1995, Eq. (B37)
            temp_mult1 = bdgfu[icl,:,iphi]*uuu[icl]*stdns[icl]/cthetv # for code readibility, multiply first the arrays with length nth
            bc1u = temp_mult1[..., np.newaxis]*gammout/cthets
            #  WWW this should be Eq. (3) of manual ver. 09.2002. Since the only variable
            #    depending on x,y,z is bidirectional gap probability, it is integrated
            #    beforehand (bdgfu). but WHY IS IT DIVIDED BY cthetv? Is it because of
            #    the definition of Gamma?

            # Calculate contribution of trunk reflectance
            # Gamma, vertical orientation (tty = 0, rfr = 0, upper hemisphere)
            gammout = TrunkRefl[icl]*(sphi - ( phi - np.pi)*cphi)/(2*np.pi)
            temp_mult2 = btr1uk[icl,:,iphi]*stdns[icl]*sthetv/cthetv
            btr1u = temp_mult2[..., np.newaxis]*gammout*tgths
            bi0u[:,iphi,:]  +=  bc1u + btr1u

            # ------ lower hemisphere
            gammd = gammt[..., np.newaxis]*LeafBranchReflLamb[icl] + gammr[..., np.newaxis]*LeafBranchTrans[icl]
            rfr = gmfres(calp2b[..., np.newaxis], RefrIndex[icl], pkhair )
            gammout = gammd + rfr*.125
            temp_mult3 = bdgfd[icl,:,iphi]*uuu[icl]*stdns[icl]/cthetv
            bc1d = temp_mult3[..., np.newaxis]*gammout/cthets
            #  Gamma,  spherical orientation, lower hemisphere
            gammout = TrunkRefl[icl]*(sphi - phi*cphi)/(2*np.pi)
            temp_mult4 = btr1dk[icl,:,iphi]*stdns[icl]*sthetv/cthetv
            btr1d = temp_mult4[..., np.newaxis]*gammout*tgths

            bi0d[:,iphi,:]  +=  bc1d + btr1d
        # end loop over tree classes
    return bi0u, bi0d
#   --------------------------- end hetk8o


def gmfres(CosAlpha, RefrIndex, pkhair ):
    """ Gamma function Fresnel reflection

    Original fortran 77 code by A. Kuusk (02.01.1991), python translation Matti Mõttus (2022)
    Computes the gamma-coeficient for Fresnel reflection at a boundary of media with
    a difference in refractive indices of RefrIndex, which depends only on the cosine
    between the incidence angle alpha relative to boundary normal. The cosine of the
    exit direction is the same as that of the incidence direction for Fresnel reflection.

    Args:
    CosAlpha: cosine of the angle of incident radiation, cos(th_incident) (f77 calp2)
    RefrIndex: refract_ind (f77 rn)
    pkhair = leaf hair index

    Returns:
    Gamma-function of reflectance? (f77 gmf)
    """
    x2  = CosAlpha*CosAlpha
    ag  = 2*x2 - 1 + RefrIndex*RefrIndex
    bg  = 1 + (ag - 2)*x2
    xy  = ag - x2
    cg  = 2*CosAlpha*np.sqrt(xy)
    sa2 = 1 - x2
    yg  = (bg + sa2*cg)*(ag + cg)
    yg  = (ag - cg)*bg/yg
    yy  = np.sqrt(sa2)/(np.pi/2)/CosAlpha*pkhair
    return np.exp(-yy)*yg
#   -------------------------- end gmfres


def pi11u(tgthx, zz12, htree, cellb, rcri):
    """Projection of the part of ellipsoid above zz12 on a horizontal surface

    translated from f77 code by A. Kuusk (16.11.2000)
    Ellipsoidi projektsioon horisontaalsele pinnale szux
    ko~rgusel zz12 suunas thx ja u"lemise osa ruumala vzui,
    tgthx = tan(thx), htree - puu ko~rgus,
       cellb - ellipsi pooltelg, rcri - krooni raadius

    Args:
    tgthx: tan(t), where t = zenith angle of the ray used for projection
    zz12: height of the surface
    htree: tree height
    cellb: vertical semiaxis of the ellipsoid (crown length/2)
    rcri: horizontal semiaxis of the ellipsoid (crown radius)

    Returns:
    szux: projected area of the upper part of the ellipsoid
    vzui: volume of the upper part of the ellipsoid
    """
    eps = 0.1e-4

    if (zz12 >= htree):
        # zz12 is above the crown
        return 0,0
    hcrown = 2.0*cellb
    hbase  = htree - hcrown

    if (zz12 >= 0 and zz12 <= hbase):
        # zz12 is below the crown
        # return total projected area and total crown volume
        szux = np.sqrt(rcri**2 + (cellb*tgthx)**2)*np.pi*rcri
        vzui = 4*np.pi*rcri**2*cellb/3
        return szux, vzui

    hcent = htree - cellb # height of the centre of the ellipsoid
    zx    = zz12 - hcent # height relative to crown center
    vzui  = np.pi*rcri**2*(hcrown/3 - zx + zx**3/3/cellb**2)
    # MATTI WWW: vzui could occasionally become slightly negative causing problems downstream
    if (vzui < eps):
        vzui=0
    # volume of crown above zz12
    #     puutepunkti kõrgus
    #     height above crown center at which the ray touches the ellipsoid
    if (tgthx < eps):
        z0 = 0
    else:
        z0 = cellb/np.sqrt(rcri/(cellb*tgthx)**2 + 1)
    if (cellb**2 - zx**2) < 0: # XXX DEBUG. this if-clause should be removed once problem is clarified
        print("XXX ERROR in pi11u, (cellb**2 - zx**2) < 0 ")
        print("XXX z12, hbase, htree, cellb, zx")
        print(z12, hbase, htree, cellb, zx)
    rz = np.sqrt(cellb**2 - zx**2)*rcri/cellb # WWW this line can throw an exception
    if ( np.abs(zx) >= z0 ):
        if (zz12 > hcent):
            # võra ülemine ots
            szux = np.pi*rz**2
        else:
            # võra alumine ots
            szux = np.sqrt(rcri**2 + (cellb*tgthx)**2)*np.pi*rcri
    else:
        # puutepunktide vahel
        xyz  = np.sqrt(rcri**2 + (cellb*tgthx)**2)
        beta = np.arccos(zx*tgthx/xyz)
        sel1 = rcri*xyz*beta - rz*zx*tgthx
        sel2 = np.pi*rz**2/2
        szux = sel1 + sel2
    return szux, vzui
#   -------------------------- end pi11u


def pi11d(tgthx, zz12, htree, cellb, rcri):
    """ projection and volume of the lower part of the ellipsoid

    subroutine pi11d, A. Kuusk 16.11.2000; python translation Matti Mõttus (2022)
    Ellipsoidi alumise osa projektsioon szdx ja ruumala vzdi
    ko~rgusel zz12 suunas thx
    tgthx = tan(thx), htree - puu ko~rgus,
    cellb - ellipsi pooltelg, rcri - krooni raadius

    Args:
    tgthx: tan(thx), thx = direction of projection
    zz12: cutoff(?) height
    htree: tree height (=crown base height + crown height)
    cellb: vertical half-axis of the ellipsoid (=1/2 crown height)
    rcri: crown (ellipsoid) radius

    Returns:
    szdx: projection area of the lower part of the ellipsoid
    vzdi: volume of the lower part of the ellipsoid
    """
    eps = 0.1e-4
    hcrown = 2*cellb
    hbase  = htree - hcrown
    #                                                trunk, tüvi
    if (zz12 >= 0) and (zz12 <= hbase):
        return 0, 0 # we are above crown, võrast kõrgemal
    if zz12 >= htree:
        vzdi = 4*np.pi*rcri**2*cellb/3
        szdx = np.sqrt(rcri**2 + (cellb*tgthx)**2)*np.pi*rcri
        return szdx, vzdi

    hcent = htree - cellb
    zx    = zz12 - hcent
    vzdi  = np.pi*rcri**2*(hcrown/3 + zx - zx**3/3/cellb**2)
    # MATTI: vzdi can become very slightly negative throwing an NaN downstream in f77
    if vzdi < eps:
        vzdi = 0
    #                                              puutepunkti kõrgus
    if tgthx < eps:
        z0 = 0
    else:
        z0 = cellb/np.sqrt(rcri/(cellb*tgthx)**2 + 1)

    rz = np.sqrt(cellb**2 - zx**2)*rcri/cellb
    if abs(zx) >= z0:
        if zz12 > hcent:
            #          võra ülemine ots
            szdx = np.sqrt(rcri**2 + (cellb*tgthx)**2)*np.pi*rcri
        else:
            #          võra ülemine ots
            szdx = np.pi*rz**2
    else:
        #                       puutepunktide vahel
        xyz  = np.sqrt(rcri**2 + (cellb*tgthx)**2)
        beta = np.arccos(zx*tgthx/xyz)
        sel1 = rcri*xyz*beta - rz*zx*tgthx
        sel3 = rcri*xyz*(np.pi - beta) + rz*zx*tgthx
        szdx = sel1 + sel3
    return szdx, vzdi
#   --------------------------  end pi11d


def pi22u (tgthx, zz12, htree, hc1i, hc2i, rcri):
    """ The volume and proj. area of the crown part above level zz12 for cone+cylinder

    Translated from the f77 code by A. Kuusk (14.11.2000). Tasandist zz12
    kõrgemale jääva võra osa ruumala ja projektsioon suunas thx, koonus + silinder
    tgthx = tan(thx), htree - puu ko~rgus, hc1i, hc2i - koonuse ja silindri ko~rgus, rcri - krooni raadius
    szui - ülemise osa projektsiooni pindala, vzui - ülemise osa ruumala

    Args:
    tgthx: tan(t), where t = zenith angle of the ray used for projection
    zz12: height of the surface
    htree: tree height
    hc1i: height of cone
    hc2i: height of cylinder
    rcri: crown radius

    Returns:
    szux: projected area of the upper part of the ellipsoid
    vzui: volume of the upper part of the ellipsoid
    """

    if (zz12 >= htree):
        # z > htree
        return 0, 0

    vcrown = rcri**2*np.pi*(hc1i/3 + hc2i)
    scon = scone(hc1i, hc1i, rcri, tgthx)
    #                                                tüvi
    if (zz12 < (htree - (hc1i + hc2i))):
        vzui = vcrown
        szux = 2*hc2i*rcri*tgthx + scon
        return szux, vzui
    #                                                silinder
    hcbase = htree - hc1i
    if (zz12 < hcbase):
        vzui  = rcri**2*np.pi*(hc1i/3 + hcbase - zz12)
        scylu = 2*rcri*tgthx*(hcbase - zz12)
        szux  = scylu + scon
        return szux, vzui
    #                                                koonus
    hcz  = htree - zz12
    rcz  = hcz*rcri/hc1i
    vzui = rcz**2*np.pi*hcz/3
    szux = scone(hcz, hcz, rcz, tgthx)

    return szux, vzui
#   --------------------------- end pi22u


def pi22d (tgthx, zz12, htree, hc1i, hc2i, rcri):
    """ Volume and projected area of the part below zz12 of the cylinder+cone crown

    subroutine pi22d   A. Kuusk   14.11.2000
    Tasandist zz12 madalamale jääva võra osa
    ruumala ja projektsioon suunas thx, koonus + silinder
    tgthx = tan(thx), htree - puu ko~rgus,
    hc1i, hc2i - koonuse ja silindri ko~rgus, rcri - krooni raadius
    szdi - alumise osa projektsiooni pindala, vzdi - alumise osa ruumala

    Args:
    tgthx: tan(thx), thx = view (projection) zenith angle
    zz12: height of the cutoff level
    htree: tree height (=crow height + crown base height)
    hc1i: length of the cone
    hc2i: length of the cylinder
    rcri: crown radius

    Returns:
    szdi: projected area of the lower part
    vzdi: volume of the lower part
    """

    hbase  = htree - (hc1i + hc2i)
    #                                         stem, tüvi
    if zz12 <= hbase:
        return 0, 0
    vcrown = rcri**2*np.pi*(hc1i/3.0 + hc2i)
    scyl   = 2.0*hc2i*rcri*tgthx
    #                                          z > htree
    if zz12 >= htree:
        vzdi   = vcrown
        scon = scone(hc1i, hc1i, rcri, tgthx)
        scrown = scyl + scon
        szdx   = scrown
        return szdx, vzdi
    # cylinder, silinder
    hcbase = htree - hc1i
    if zz12 < hcbase:
        vzdi  = rcri**2*np.pi*(zz12 - hbase)
        scyld = 2.0*rcri*tgthx*(zz12 - hbase)
        szdx  = scyld + rcri**2*np.pi
        return szdx, vzdi
    #  cone, koonus
    hcz  = htree - zz12
    rcz  = hcz*rcri/hc1i
    vzui = rcz**2*np.pi*hcz/3.0
    vzdi = vcrown - vzui
    scdx = scone(hc1i, hcz, rcri, tgthx)
    szdx = scdx + scyl
    return szdx, vzdi
#   -------------------------- end pi11d


def scone(hcone, hcz, rcone, tgthx ):
    """ (Tüvi)koonuse projektsioon
    Projection of a frustum

    translated from f77 code by A. Kuusk (14.11.2000)
    hcone - koonuse kõrgus, hcz - tüvikoonuse kõrgus, rcone - koonuse
    aluse raadius, tgthx - vaatenurga tangens

    Args:
    hcone: height of the cone forming the frustum
    hcz: height of the frustum
    rcone: radius of the base
    tgthx: tan(t), where t = zenith angle of the ray used for projection

    Returns:
    scz: projection area
    """

    eps = 0.1e-4
    rhc   = rcone/hcone/np.max( (eps, tgthx) )
#                                              thx > thc
    if (rhc < 1):
        hc4  = np.abs(hcone - hcz)
        beta = np.arccos(rcone/(hcone*tgthx))
#                                        tüvikoonus
        if (hc4 > eps):
            rcz = rcone*hc4/hcone
            sc3 = rcz**2*beta
            sc4 = np.sqrt((hc4*tgthx)**2 - rcz**2)*rcz
        else:
            sc3 = 0
            sc4 = 0

        sc1 = np.sqrt((hcone*tgthx)**2 - rcone**2)*rcone
        sc2 = rcone**2*(np.pi - beta)
        scz = sc2 + sc1 - sc4 + sc3
    else:
        scz  = np.pi*rcone**2
    return scz
#   -------------------------- end scone


# =============================================================
#  functions for twostream computations


def twostr(thets, thetv, sqratio, efflai, tlty, cf, rteff, tteff, rdsgrou,
     rsdgrou, rddgrou, t0, tv ):
    """ Two-stream calculations for calculating diffuse fluxes in FRT

    Part of original frt by A. Kuusk, modified and translated to python by Matti Mõttus (2022)

    Args:
    thets: solar zenith, rad
    thetv: view zenith, rad
    efflai: Canopy LAI. In FRT, an effective (optical) LAI value of leaves + branches is used
    sqratio: radio of direct to total irradiance
    rteff: mean Lambertian component of the spectral reflectance of leaves+branches
    tteff: mean (diffuse) transmitance of leaves+branches
    tlty: stem area index, separated from leaves + branches as it has a different orientation
    rtrnk: mean trunk reflectance
    cf: correction factor for 1st and higher order scattering
    leaf hair optical index (vaguely defined, set to unity), input parameter for FRT
    rdsgrou: ground hemispherical-directional reflectance (diffuse incidence)
    rsdgrou: ground directional-hemispherical reflectance
    rddgrou: ground bi-hemispherical reflectance (diffuse incidence) = albedo
    t0: gaps_s(1) from hetk8s: transmittance in solar direction
    tv: gaps_s(1) from hetk8s: transmittance in view direction

    Returns:
    Rhd_hi_c: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-direcitpnal),
        contribution by canopy, includes all-order diffuse-sky reflectance component
    Rhd_hi_g: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-direcitpnal),
        contribution by ground (forest floor), includes all-order diffuse-sky reflectance component
    Thd_hi: higher-order component (excl. 1st-order direct) of canopy+soil hemispherical-directional transmittance
        includes all-order diffuse-sky reflectance component
    the higher-order components (*_hi) include all orders of reflectance of diffuse-sky radiation
    """

    # Notation: [reflectance,transmittance][hemispherical,directional,isotropic][directional,hemispherical]_[higherorder,total][Blacksurface]_[canopy,ground]
    # tdi_t: total directional-isotropic layer scattering factor (tsd in Kuusk 2001) FOR BLACK SOIL
    # rid_t: total isotropic-directional layer scattering factor (rsd Kuusk 2001) FOR BLACK SOIL
    # rii_t: total isotropic-isotropic layer scattering factor (rdd Kuusk 2001) FOR BLACK SOIL
    # tid_t: total isotropic-directional layer transmittance factor (tsd Kuusk 2001) FOR BLACK SOIL
    # tii_t: total isotropic-isotropic layer transmittance factor (tdd Kuusk 2001) FOR BLACK SOIL

    gsf = 0.5 # spherical G
    skyl  = 1 - sqratio

    # Compute layer scattering factors (i.e., black ground)
    utot = efflai + tlty # lai + trunk area index
    gcyl = 2*np.sin(thets)/np.pi
    ggef = (efflai*gsf + tlty*gcyl)/utot
    #  ggj  - the J-integral in the SAIL model multiplied by ul_z / 2
    ggj  = efflai/6/utot
    rii_t, tii_t, rid_t, tid_t, rdi_t, tdi_t, rdd_hi = \
        layer( thets, thetv, ggef, ggef, ggj, utot, rteff, tteff )
    # NOTE: layer does not output tsot, sun->observer diffuse transmittance.
    # NOTE: rdi_t not used below

    # Compute R,T which includes the effect of reflective ground
    #   correction factor is only applied to directional->X quantities
    fac  = 1/(1 - rii_t*rddgrou)

    Rhd_hi_c   = sqratio*cf*rdd_hi + skyl*rid_t \
        + (sqratio*t0*rsdgrou + sqratio*cf*tdi_t*rddgrou \
        + skyl*tii_t*rddgrou)*tid_t*fac        # Kuusk 2001,  (9)

    # Rhd_hi_g = (sqratio*tdi_t + skyl*tii_t)*rdsgrou*tv*fac
    #  Kuusk 2001, (10) with added fac -- probably en error in the original eqn
    # Modified in 2022 Rhd_hi_g to include the contribution by direct
    # transmittance + multiple scattering
    Rhd_hi_g=(sqratio*cf*tdi_t + skyl*tii_t + sqratio*t0*rsdgrou*rii_t)*rdsgrou*tv*fac

    # Compute higher-order down in the direction of observation,
    # AKA higher-order directional transmittance, Thd_hi
    # bgrff can output Thd_hi for a reflective ground, but its output has not been corrected
    #  with cf  and the subroutine is difficult to modify.
    # Hence, use bgrdd to compute tdd_hi (a.k.a. tsot in Kuusk 2001), the layer scattering
    # coefficient for higher-order scattering, and add the effect of reflecting surface later.
    #  To compute tdd_hi, set skyl and ground reflectance to zero
    tdd_hi = bgrdd(thets, thetv, 0, efflai, tlty, rteff, tteff, 0, 0)

    # compose Thd_hi from the layer scattering factors
    Thd_hi = sqratio*cf*tdd_hi + skyl*tid_t + \
        ( sqratio*(t0*rsdgrou+tdi_t*cf*rddgrou) +  \
        skyl*tii_t*rddgrou ) *rid_t / (1-rii_t*rddgrou)
    return Rhd_hi_c, Rhd_hi_g, Thd_hi
#   --------------------------- end twostr


def layer(thets, thetv, ggs, ggv, ggj, ul, rrl, ttl):
    """ the scattering operators of a leaf layer

    Based on the fortran 77 code by A. Kuusk (13.01.2000), Diffuse reflection and
    transmission of a layer from the 2-stream model for elliptical (?) leaf orientation.
    Calculates the scattering operators of a leaf layer (Table 1; Kuusk 2001)
    Code is similar to that given in Kuusk 1995, A computer-efficient...
        Comp & Geosci 22, 149-163, Appendix B
      and Kuusk 1995, A fast invertible..., RSE 51, 342-350

    Other outputs than rso include most probably all scattering orders, but this is contradictory in (Kuusk, 2001)
      Kuusk (2001) Eq. (9) for rho_d^plants includes SQr_so, although it should not include by definition
      first-order scattering from the direct beam. Exclusion of first-order scattering from r_so is not
      explicitly stated anywhere, but is evident from Eq. (9).

    Args:
       thets - the Sun zenith, thetv - view angle
       ggs, ggv - G-function at thets and thetv
       ggj  - the J-integral divided by 2 (Eq. 6)
       ul - LAI
       rrl, ttl - leaf reflection and transmission

    Returns: rdd, tdd, rdo, tdo, rsd, tsd, rso
       Notation: [r=reflectance,t=transmittance][s=sun,d=diffuse][o=observer,d=diffuse]
       rso - higher-order reflectance of the layer in view direction
             for direct incoming radiation
    """

    rtp   = (rrl + ttl)/2.0
    rtm   = (rrl - ttl)/2.0
    bf    = rtm*ggj

    #  Eq. numbers in Kuusk 1996, A computer-efficient, Comp & Geosci 22, 149-163, App. B
    vks   = ggs/np.cos(thets)        # k_1 from B30, excl u_L
    vkv   = ggv/np.cos(thetv)        # k_2 from B30, excl u_L
    vsig  = rtp + bf                       # B23, excl u_L
    vatt  = 1.0 - rtp + bf                # B22, excl u_L
    vssf  = vks*rtp - bf                   # B24
    vssb  = vks*rtp + bf                   # B25
    vmm   = np.sqrt(vatt**2 - vsig**2)        # B21
    vvv   = vkv*rtp + bf
    vuu   = vkv*rtp - bf
    vh1   = (vatt + vmm)/vsig              # B20
    vh2   = 1.0/vh1
    vcc   = (vssf*vsig - vssb*(vks - vatt))/(vmm**2 - vks**2) # B16
    vdd   = (vssb*vsig + vssf*(vks + vatt))/(vmm**2 - vks**2)

    exmu1 = np.exp(vmm*ul) # ???
    exmu2 = 1.0/exmu1
    exku1 = np.exp(-vks*ul)
    delt  = vh1*exmu1 - vh2*exmu2

    rdd   = (exmu1 - exmu2)/delt
    tdd   = (vh1 - vh2)/delt

    delta = vh2*vcc*exku1 - vdd*exmu1
    aa1   = delta/delt
    deltb = vdd*exmu2 - vh1*vcc*exku1
    bb1   = deltb/delt
    tsd   = vh1*aa1*exmu2 + vh2*bb1*exmu1 + vdd*exku1
    rsd   = aa1 + bb1 + vcc
    tdo   = ul*((vuu*vh1 + vvv)*eint((vkv - vmm)*ul) -
            (vuu*vh2 + vvv)*eint((vkv + vmm)*ul))/delt
    rdo   = ul*(exmu1*(vvv*vh1 + vuu)*eint((vkv + vmm)*ul) -
            exmu2*(vvv*vh2 + vuu)*eint((vkv - vmm)*ul))/delt
    rso   = ul*((vvv*vdd + vuu*vcc)*eint((vkv + vks)*ul) +
            (vvv*vh1 + vuu)*aa1*eint((vkv + vmm)*ul) +
            (vvv*vh2 + vuu)*bb1*eint((vkv - vmm)*ul))
    return rdd, tdd, rdo, tdo, rsd, tsd, rso
#   -------------------------- end layer


def bgrdd(thets, thetv, skyl, efflai, tlty, rteff, tteff, rddgrou, rsdgrou):
    """ Twosteam canopy radiation model, diffuse (higher-order) fluxes down

    Based on Fortran77 subroutine bgrdd, Difuussed vood 2-voo lähenduses
    A. Kuusk  14.10.2000, 12.09.2002; translation to python Matti Mõttus 2022

    Args:
    thets, thetv, skyl,
    efflai: efffective (optical) LAI+BAI
    tlty: LAI of trunks
    tlty, rteff, tteff: optical properties
    rddgrou: bi-hemispherical ground reflectance (diffuse incidence)
    rsdgrou: directional-hemispherical ground reflectance

    Returns:
    bddif: diffuse fluxes down
    """

    gsf = 0.50 # G = 0.5, spherical leaf orientation
    ul    = efflai + tlty
    gcyl  = 2.0*np.sin(thets)/np.pi
    ggef  = (efflai*gsf + tlty*gcyl)/ul
    ggs   = ggef
    ggv   = ggef
    ggj   = efflai/6.0/ul

    rtp  = (rteff + tteff)/2.0
    rtm  = (rteff - tteff)/2.0
    bf   = rtm*ggj
    #     cthetv = cos(thetv)
                                # v??? - Verhoefi tähistused
    vks   = ggs/np.cos(thets)
    vkv   = ggv/np.cos(thetv)
    vsig  = rtp  + bf
    vatt  = 1.0 - rtp + bf
    vssf  = vks*rtp - bf
    vssb  = vks*rtp + bf
    vmm   = np.sqrt(vatt**2 - vsig**2)
    vvv   = vkv*rtp + bf
    vuu   = vkv*rtp - bf
    vh1   = (vatt + vmm)/vsig
    vh2   = 1.0/vh1
    vcc   = (vssf*vsig - vssb*(vks - vatt))/(vmm**2 - vks**2)
    vdd   = (vssb*vsig + vssf*(vks + vatt))/(vmm**2 - vks**2)

    exmu1 = np.exp(vmm*ul)
    exmu2 = 1.0/exmu1
    exku1 = np.exp(-vks*ul)

    delt  = exmu2*(vh2 - rddgrou) - (vh1 - rddgrou)*exmu1
    delta = (1.0 - skyl)*exku1*(vdd*rddgrou + rsdgrou - vcc)*vh2 - \
            (1.0 - rddgrou*vh2)*(skyl - vdd*(1.0 - skyl))*exmu1
    deltb = exmu2*(1.0 - vh1*rddgrou)*(skyl - vdd*(1.0 - skyl)) - \
            vh1*(1.0 - skyl)*exku1*(vdd*rddgrou + rsdgrou - vcc)
    vaa   = delta/delt
    vbb   = deltb/delt

    bddif = vaa*(vuu*vh1 + vvv)*ul*exmu2*eint(ul*(vkv - vmm)) + \
            vbb*(vuu*vh2 + vvv)*ul*eint(ul*(vkv + vmm))*exmu1 + \
            (vuu*vdd + vvv*vcc)*ul*(1.0 - skyl)*exku1* \
            eint(ul*(vkv - vks))
    return bddif
#   -------------------------- end bgrdd
