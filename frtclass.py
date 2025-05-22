#!/usr/bin/env python3
# the FRT class: forest reflectance and transmittance. Before executing, check the paths at the end of the file,
#  especially frt_datadir. Executing this code creates the object G of type frt_model, configures it,
#  and computes it's reflectance and transmittance. The results are in G.R and G.T with more details
#  -- such as scattering by different orders -- in other variables of G.
#  to compute flux reflectance, run G.flux_reflectance(). See the relevant functions for their outputs.
#
import numpy as np
import os
import sys
from scipy.optimize import brentq
from copy import deepcopy
import json

from frtfunctions import *
from frtfunctions_py import *
# to use the fortran77 bindings and libraries, use this instead of the above
# from frtfunctions_f77 import *
# NOTE! to use fortran modules, run first "make frt" in the source code directory
# NEXT, the following modules need to be compiled with f2py (see their individual .f files for details)
#   spooi, enel3, bck3

# to use f-art input files, run also this (or similar, depending on system config)
#  in the frt source code folder:
# f2py.exe -c --compiler=mingw32 -m xd_cfm xd_cfm.f


class frt_model:

    def __init__( self, frt_srcdir=None ):
        # frt_srcdir is needed to import the fortran modules
        self.frtconf = {}
        self.frtconf_isread = False
        self.configuration_applied = False
        self.frt_configured = False
        self.frt_srcdir = frt_srcdir # frt_srcdir is None: python-only FRT.
        if frt_srcdir is not None:
            if frt_srcdir not in sys.path:
                sys.path.append(frt_srcdir)
        self.frt_datadir = "."
        # the required key lists below are not complete and do not guarantee a minimum set
        self.NeededConfKeys = ['TreeClasses', 'thetv', 'phiv', 'thets', 'wl']
        self.NeededClassKeys = ['l_elli', 'StandDensity', 'TreeHeight', 'CrownLength1',
            'CrownRadius', 'DBH', 'DryLeafWeight', 'SLW', 'BAILAI', 'TreeClumping', 'ShootClumping',
            'ShootLength' ]
        self.correctFGI = True # whether to correct for the unrealistic case of canopy cover  > crown cover
        self.correctcrownlength = True # whether to correct for unrealisticly long crowns

    def load_conf( self, frtconf, frt_datadir=None ):
        """ Load the external dict frtconf into frt and do some basic checks, load spectrum files
        optionally, path to data files may be given
        """
        if frt_datadir is not None:
            self.frt_datadir = frt_datadir
        for key in self.NeededConfKeys:
            if key not in frtconf.keys():
                print("frt_model.load_conf(): Required key "+key+" missing in frtconf. Aborting")
                return
        for i,tc in enumerate(frtconf['TreeClasses']):
            for key in self.NeededClassKeys:
                if key not in tc.keys():
                    print("frt_model.load_conf(): Required key "+key+" missing in TreeClass "
                    + str(i) + ". Aborting")
                    return
        # check existence of spectral data
        self.frtconf = deepcopy(frtconf) # just in case :), we will modify frtconf while loading
        wl = self.frtconf['wl']
        SpectralKeys = [ 'rDDground', 'rDSground', 'rSDground', 'rSSground', 'WaxRefrInd', 'SQratio']
        for key in SpectralKeys:
            if key not in self.frtconf.keys(): # check for data file name
                if key+'File' not in self.frtconf.keys():
                    print("frt_model.load_conf(): Required "+key+" or "+key
                        +"File not found in frtconf. Aborting." )
                    return
                else: # read data from file and interpolate
                    fname = os.path.join( self.frt_datadir, self.frtconf[ key+'File' ] )
                    print("reading "+fname)
                    Q = np.genfromtxt( fname )
                    wl_Q = Q[:,0]
                    v_Q = Q[:,1]
                    self.frtconf[key] = np.interp( wl, wl_Q, v_Q )
                #
            #
        # Next, do the same for all tree classes
        for i,tc in enumerate(self.frtconf['TreeClasses']):
            SpectralKeys = [ 'LeafRefl', 'LeafTrans', 'BranchRefl', 'TrunkRefl' ]
            for key in SpectralKeys:
                if key not in tc.keys(): # check for data file name
                    if key+'File' not in tc.keys():
                        print("frt_model.load_conf(): Required "+key+" or "+key
                            +"File for TreeClass " +int(i)+" not found in frtconf. Aborting." )
                        return
                    else: # read data from file and interpolate
                        fname = os.path.join( self.frt_datadir,tc[ key+'File' ] )
                        print("reading "+fname)
                        Q = np.genfromtxt( fname )
                        wl_Q = Q[:,0]
                        v_Q = Q[:,1]
                        tc[key] = np.interp( wl, wl_Q, v_Q )
                        if key=='LeafRefl':
                            # Check if Q has more than two columns
                            if Q.shape[1] > 2:
                                # leaf transmittance in column 3
                                v_Q = Q[:,2]
                                tc["LeafTrans"] = np.interp( wl, wl_Q, v_Q )
                                SpectralKeys.remove('LeafTrans')
                                if Q.shape[1] > 3:
                                    # abaxial reflectance in column 4
                                    v_Q = Q[:,3]
                                    tc['LeafRefl2'] = np.interp( wl, wl_Q, v_Q )
                                #
                            #
                        #
                    #
                #
            # Finally, the fully optional (and not documented?) leaf abaxial reflectance as a separate file
            if 'LeafRefl2File' in tc.keys():
                fname = os.path.join( self.frt_datadir,tc['LeafRefl2File'] )
                print("reading "+fname)
                Q = np.genfromtxt( fname )
                wl_Q = Q[:,0]
                v_Q = Q[:,1]
                tc['LeafRefl2'] = np.interp( wl, wl_Q, v_Q )
            if 'ScaleNeedle' in tc.keys():
                if tc['ScaleNeedle']:
                    # we need to scale the reflectance from needle -> shoot
                    #  using basic p-theory, p computed from ShootClumping
                    p = 1 - tc['ShootClumping']
                    # calculate leaf albedo
                    if 'LeafRefl2' in tc.keys():
                        R = np.array( tc['LeafRefl'] + tc['LeafRefl2'] )/2
                    else:
                        R = np.array( tc['LeafRefl'] )
                    w = R + np.array( tc['LeafTrans'] )
                    # preserve R/T ratio. Currently, no better option.
                    tc['LeafRefl'] = np.array( tc['LeafRefl'] )*(1-p)/(1-p*w)
                    tc['LeafTrans'] = np.array( tc['LeafTrans'] )*(1-p)/(1-p*w)
                    if 'LeafRefl2' in tc.keys():
                        tc['LeafRefl2'] = np.array( tc['LeafRefl2'] )*(1-p)/(1-p*w)
                    print("Scaled leaf R & T using p=", str(p)," for tree class ",
                        str(i), end=" ")
                    if 'Description' in tc.keys():
                        print( tc['Description'] )
                    else:
                        print("")
                    #
                #
            #
        #

        self.frtconf_isread = True
        self.configuration_applied = False # the data from the dict needs to be copied into class variables
        self.frt_configured = False # new configutration read, indicate that preparatory computations are not done

    def read_conf( self, configfilename, frt_datadir=None ):
        # read a FRT configuration file into the frtconf dictionary
        #  currently, there is no python way to read the input file, only via the fortran module
        if frt_datadir is not None:
            self.frt_datadir = frt_datadir
        # import python library
        import xd_cfm
        # import also the helpers to convert data from fortran
        from frt_wrapper_functions import chararray2strarray
        # xd_cfm needs to be run in the  folder where the data are
        os.chdir(self.frt_datadir)

        q=xd_cfm.xd_cfm( configfilename )
        # subroutine xd_cfm( fname, XC, IC, SPCIN, DESC, lerr)
        # q: ( XC, IC, SPCIN, DESC, lerr)
        XC = q[0]
        IC = q[1]
        SPCIN = q[2]
        DESC = q[3]
        lerr = q[4]

        if q[4] != 0:
            #error has occurred and error message has been  displayed.
            print("WARNING! xd_cfm detected an error in the file, read_conf() ABORTING.")
            print("xd_cfm  error flag ", lerr)
            return

        # convert  DESC into readable form. NOTE: this may depend on the specific implementation of f2py
        # DESC = chararray2strarray(DESC) # OBSOLETE?
        DESC = [ x.decode('utf8').strip() for x in DESC ] # decode and strip to get an array of proper strings
        self.frtconf["Description"] = DESC[0]
        nclmax = xd_cfm.frtpar.fnclmax[()] # extract from zero-dimensional ndarray

        # Move configuration data to the frtconf dictionary, start with integers
        #   IC elements === NOTE: BELOW ARE THE FORTRAN NUMBERS STARTING FROM 1. Subtract 1 for python
        #  1: no. of canopy classes (old n_cl)
        #  2: correction wavelength index
        #  3: no. of crown layers in numerical integration (old n_zs)
        #  4: logical flag l_opt    : whether to (re)compute element mean optical properties
        #  5: logical flag l_struc : whether to (re)compute canopy overall structural characteristics
        #  6: logical flag l_dirv   : whether to (re)compute everything related to view direction (gaps etc.)
        #  7: logical flag l_dirs   : whether to (re)compute everything related to solar direction (gaps etc.)
        #  8: logical flag l_grnd   : whether to (re)load forest floor spectra
        #  9: logical flag l_refl   : whether to compute reflectance factors; if false,
        #         the flags are not input parameters and will be ignored here
        #  10 & 11: start and end wavelength indices, respectively. Indices start from 1 (fortran77 convention)
        #  11+i to 11+nclmax: are crown shapes ellipsoids? if >0, then yes (old l_elli)
        #  12+nclmax: ijob (used to select computation type -- not all settings can be used in python)
        #  13+nclmax: nquad_t, quadrature knots over theta (polar, zenith angle)
        #  14+nclmax: nquad_p, quadrature knots over phi (azimuth)

        # Save forest parameters in the dictionary frtconf
        #   model parameters as class parameters
        # some parameters may be both
        n_classes = IC[0]
        self.frtconf["TreeClasses"] = []
        for i in range( n_classes ):
            self.frtconf["TreeClasses"].append( { "Description": DESC[i+1] } )
        # read the correction wavelength index. The subroutine always returns a positive integer.
        #   For no correction, fortran code sets it as "IC(2) = IC(11)+100"
        # in python, integer indices are gien as negative values
        self.frtconf["wlcorr"] = -IC[1]
        self.frtconf["nzs"] = IC[2]
        # the flags below (IC[4:9]) are irrelevant and are not read from the input file
        self.frtconf["i1"] = IC[9] # retain indexing from 1, switch to python style later
        self.frtconf["i2"] = IC[10] # retain indexing from 1, switch to python style later
        #  all wavelength-oriented arrays will be used in the [i1,i2] interval. While the original data
        #   are preserved in the frtconf dictionary, the actual variables (e.g., self.wl) will be cropped to [i1,i2]
        for tc_dict in self.frtconf["TreeClasses"]:
            tc_dict["l_elli"] = IC[11+i] == 1
        self.frtconf["ijob"] = IC[11+nclmax]
        self.frtconf["nquad_t"] = IC[12+nclmax]
        self.frtconf["nquad_p"] = IC[13+nclmax]

        #  XC elements. Column refers to 2nd index
        #      === NOTE: BELOW ARE THE FORTRAN NUMBERS STARTING FROM 1. Subtract 1 for python
        # StandDensity: stand density [m^-2], XC column 1; f77 stdns
        # TreeHeight: tree height [m], XC column 2;  f77 htr
        # CrownLength1: crown length, ellipse or conical part [m], XC column 3; f77 hc1
        # CrownLength2: crown length, cylinder [m], XC column 4; f77 hc2
        # CrownRadius: crown radius [m], XC column 5; f77 rcr
        # DBH: trunk diameter [cm], XC column 6. When read into self.DBH, converted to [m]; f77 dbh
        # DryLeafWeight: dry leaf weight [kg/tree], XC column 7; f77 rmass
        # SLW: specific leaf weight, SLW [gm^-2], XC column 8; f77 slwcl
        # BAI/LAI ratio: XC column 9, NOT stored in a variable
        # TreeClumping: clumping index caused by tree distribution, cB in eq. (11) in Nilson 1999, XC column 10; f77 clmpst
        # ShootClumping: shoot shading coefficient (ssc), XC column 11; f77  clmpsh
        # crncl: wax refractive index ratio, XC column 12; f77  crncl
        # ShootLength: hot spot parameter, shoot length (m), XC column 13; f77 shl
        # positions 14-20 unused
        # thetv: view zenith angle [rad], XC(1,21)
        # phiv: view azimuth angle [rad], XC(1,22)
        # thets: solar zenith angle [rad], XC(1,23)
        # pkhair: leaf hair optical index (vaguely defined, set to unity), XC(1,24); f77 pkhair [was p_khair]
        # age: stand age (not used in computations) XC(1,25)
        # dthetv: thetv increment [rad] (not used in current computations) XC( 1,26);
        for i,tc_dict in enumerate(self.frtconf["TreeClasses"]):
            tc_dict["StandDensity"] = XC[i,0]
            tc_dict["TreeHeight"] = XC[i,1]
            tc_dict["CrownLength1"] = XC[i,2]
            tc_dict["CrownLength2"] = XC[i,3]
            tc_dict["CrownRadius"] = XC[i,4]
            tc_dict["DBH"] = XC[i,5] # note: in tc_dict, dbh is in cm.
            tc_dict["DryLeafWeight"] = XC[i,6]
            tc_dict["SLW"] = XC[i,7]
            tc_dict["BAILAI"] = XC[i,8]
            tc_dict["TreeClumping"] = XC[i,9] # clumping index caused by tree distribution
            #  =cB in eq. (11) in Nilson 1999
            tc_dict["ShootClumping"] = XC[i,10]
            tc_dict["WaxCorrectionFactor"] = XC[i,11] # 'crncl' in f77
            tc_dict["ShootLength"] = XC[i,12]
        self.frtconf["thetv"] = XC[0,20]
        self.frtconf["phiv"] = XC[0,21]
        self.frtconf["thets"] = XC[0,22]
        self.frtconf["pkhair"] = XC[0,23]
        self.frtconf["Age"] = XC[0,24]
        self.load_optional_confparameter( "age", 0 )
        self.frtconf["dthetv"] = XC[0,25]

        # ------  SPCIN elements (columns). Column refers to 2nd index
        #      === NOTE: BELOW ARE THE FORTRAN NUMBERS STARTING FROM 1. Subtract 1 for python
        # 1  wavelength (nm)
        # 2  hemispherical-hemispherical reflectance of forest floor (=albedo) (rddgrou)
        # 3  directional-hemispherical reflectance of forest floor (rsdgrou)
        # 4  hemispherical-directional reflectance of forest floor (rdsgrou)
        # 5  BRDF of forest floor (directional-directional reflectance) (rsogrou)
        # 6  S/Q ratio (direct/total irradiance) at TOC (s_qarr)
        # 7  leaf wax refractive index
        # 8          to 8+nclmax-1    leaf reflectance (p_rarr) for each individual tree class
        # 8+nclmax   to 8+2*nclmax-1  leaf transmittance (p_tarr)
        # 8+2*nclmax to 8+3*nclmax-1  leaf adaxial reflectance [NOT USED]
        # 8+3*nclmax to 8+4*nclmax-1  branch reflectance (b_rrarr)
        # 8+4*nclmax to 8+5*nclmax-1  trunk reflectance (t_rrarr)
        # 8+5*nclmax   correction factor for diffuse fluxes (to correct for the error in estimating 1st order scattering) (c_fact)
        #    c_fact is not input, so it's ignored here
        self.frtconf["wl"] = SPCIN[:,0]
        self.frtconf["rDDground"] = SPCIN[:,1] # Dif->Dif
        self.frtconf["rSDground"] = SPCIN[:,2] # Sun->Dif
        self.frtconf["rDSground"] = SPCIN[:,3] # Dif->Sensor
        self.frtconf["rSSground"] = SPCIN[:,4] # Sun->Sensor -- was "so" (Sun->Observer in f77 original code)
        self.frtconf["SQratio"] = SPCIN[:,5]
        self.frtconf["WaxRefrInd"] = SPCIN[:,6] # was 'rind' in f77
        for i,tc_dict in enumerate(self.frtconf["TreeClasses"]):
            tc_dict["LeafRefl"] = SPCIN[:,7+i]
            tc_dict["LeafTrans"] = SPCIN[:,7+nclmax+i]
            tc_dict["LeafRefl2"] = SPCIN[:,7+2*nclmax+i]
            tc_dict["BranchRefl"] = SPCIN[:,7+3*nclmax+i]
            tc_dict["TrunkRefl"] = SPCIN[:,7+4*nclmax+i]

        self.frtconf_isread = True
        self.configuration_applied = False # the data from the dict needs to be copied into class variables
        self.frt_configured = False # new configutration read, indicate that preparatory computations are not done
        return

    def load_optional_confparameter( self, name, defaultvalue ):
        """ Fill in the optional values based on confdata dict and a default value
        """
        if name in self.frtconf.keys():
            if name is not None:
                return self.frtconf[name]
        return defaultvalue

    def conf2json( self, filename=None, indent=4 ):
        """ dump the configuration dict as json

        Args:
        If filename is given, json is stored there
        indent: json parameter

        Returns:
        string: json dump of self.frtconf
        """
        # the nontrivial task is finding all the types json is not digesting and
        #    (numpy arrays, intc, bool_, ...) and converting them.
        # I know there are json encoders, but simpler is better.

        if self.frtconf_isread:
            fc = deepcopy( self.frtconf ) # make a copy for local editing
            i1 = fc['i1']
            i2 = fc['i2']
            for i in fc.keys():
                if isinstance( fc[i], np.ndarray):
                    fc[i] = fc[i][i1:i2+1].tolist()
                if isinstance( fc[i], np.intc ):
                    fc[i] = int( fc[i] )
            # next, go into treeclasses
            for tc in fc['TreeClasses']:
                for i in tc.keys():
                    if isinstance( tc[i], np.ndarray):
                        tc[i] = tc[i][i1:i2+1].tolist()
                    if isinstance( tc[i], np.bool_ ):
                        tc[i] = bool( tc[i] )
                    if isinstance( tc[i], np.intc ):
                        tc[i] = int( tc[i] )
            if filename is not None:
                fn = os.path.join( self.frt_datadir, filename)
                of = open( fn, "w" )
                json.dump( fc, of, indent=indent )
                of.close()
                print("saved configuration to file "+fn)
            return json.dumps(fc, indent = indent )
        else:
            print("FRT config not presebt, nothing to convert to json!")
            return None

    def conf2pickle( self, picklefile, absolutepath=False ):
        """ save the frt configuration dict in a binary pickle file.
        Default location of the file is in frt_datadir

        Args:
        picklefile: name of the file for pickled configuration
        absolutepath: if True, do not save in frt_datadir, assume picklefile is absolute
        """
        if self.frtconf_isread:
            import pickle
            if not absolutepath:
                picklefile = os.path.join( self.frt_datadir, picklefile )
            f = open(picklefile,'wb')
            pickle.dump(self.frtconf,f)
            f.close()
            print( "Saved configuration to "+picklefile )
        else:
            print("FRT config not present, nothing to pickle!")


    def load_pickle( self, picklefile, frt_datadir=None ):
        """ load the frt configuration dict from a binary pickle file.
        NOTE: the pickle file should be the one which has been already loaded by frt,
        i.e., with interpolated data and no references to external files! No external
        files will be used for optical properties etc.

        Args:
        picklefile: name of the file for pickled configuration
        frt_datadir (optional): if given, set frt_datadir -- and look here for picklefile, too
        """
        import pickle

        if frt_datadir is not None:
            self.frt_datadir = frt_datadir
        # first, search for picklefile in the datadir (i.e., assume it's not absolute name)
        if os.path.isfile( os.path.join(self.frt_datadir,picklefile) ):
            picklefile = os.path.join(self.frt_datadir,picklefile)
        self.frtconf = pickle.load( open( picklefile, "rb") )

        self.frtconf_isread = True
        self.configuration_applied = False # the data from the dict needs to be copied into class variables
        self.frt_configured = False # new configutration read, indicate that preparatory computations are not done

    def apply_config( self ):
        """copy data from conf dictionary into object
        Call this after config data in a dictionary has been loaded
        """
        if not self.frtconf_isread:
            print("Cannot apply configuration as one needs to be read/loaded first")
            return
        #
        # Compute knots and weights of Gauss-Legendre quadrature over the hemisphere
        #    nquad_t and nquad_p are the number of knots
        #    and xquad_t() and xquad_p() for the knots themselves
        # Azimuth: the model is symmetric wrt principal plane, divide nodes equally in [0,pi]
        #   if nquad_p == 1, this leads to the cross plane
        self.Description = self.load_optional_confparameter( "Description", "FRT configuration" )
        self.ncl = len( self.frtconf["TreeClasses"] )
        self.wl = np.array( self.frtconf["wl"] ).astype(float)
        self.wlcorr = self.load_optional_confparameter( "wlcorr", 0 )
        self.nzs = self.load_optional_confparameter("nzs",4)
        self.nquad_p = self.load_optional_confparameter( "nquad_p", 7 )
        self.nquad_t = self.load_optional_confparameter( "nquad_t", 7 )
        self.correctFGI = self.load_optional_confparameter( "correctFGI", True )
        self.thetv = self.frtconf["thetv"]
        self.phiv = self.frtconf["phiv"]
        self.thets = self.frtconf["thets"]
        self.dthetv = self.load_optional_confparameter( "dthetv", 5*np.pi/180 )
        self.pkhair = self.load_optional_confparameter( "pkhair", 1 )
        self.i1 = 0 # the staring index in wavelength
        if 'i1' in self.frtconf.keys():
            if self.frtconf["i1"] is not None:
                # switch to python-style indexing, starting from zero
                self.i1 = self.frtconf["i1"]-1
        self.i2 = len( self.wl )
        if 'i2' in self.frtconf.keys():
            if self.frtconf["i2"] is not None:
                # switch to python-style indexing, starting from zero
                self.i2 = self.frtconf["i2"]-1

        # store tree class parameters in arrays for computations
        self.l_elli = self.getarray('l_elli')
        self.StandDensity = self.getarray('StandDensity').astype(float)
        self.SLW = self.getarray('SLW').astype(float)
        self.DryLeafWeight= self.getarray('DryLeafWeight').astype(float)
        self.LAI = self.StandDensity*self.DryLeafWeight*1000/self.SLW # in f77: rlai
        bailai = self.getarray('BAILAI').astype(float)
        self.BAI = self.LAI*bailai
        self.ShootLength = self.getarray('ShootLength').astype(float)
        self.TreeHeight = self.getarray('TreeHeight').astype(float)
        self.CrownLength1 = self.getarray('CrownLength1').astype(float)
        self.CrownLength2 = self.getarray('CrownLength2').astype(float)
        if self.correctcrownlength:
            # Crown length should not exceed tree height. Consider the latter a
            #   more reliable parameter and adjust crown length
            for i_cl, l_cl in enumerate(self.l_elli):
                # No error is created for efficiency
                if l_cl:
                    self.CrownLength1[i_cl] = min(self.CrownLength1[i_cl],self.TreeHeight[i_cl])
                else:
                    # shrink proportionally both crown parts, cone and cylinder
                    CL = self.CrownLength1[i_cl] + self.CrownLength2[i_cl]
                    if CL > self.TreeHeight[i_cl]:
                        cf = self.TreeHeight[i_cl]/CL
                        self.CrownLength1[i_cl] *= cf
                        self.CrownLength2[i_cl] *= cf

        self.CrownRadius = self.getarray('CrownRadius').astype(float)
        self.DBH = self.getarray('DBH').astype(float)/100 # NB! convert from cm to m
        self.TreeClumping = self.getarray('TreeClumping').astype(float)
        self.ShootClumping = self.getarray('ShootClumping').astype(float)

        self.xquad_p = [ np.pi/self.nquad_p*(0.5 + i) for i in range(self.nquad_p) ]
        # The quadrature weights combine both theta and phi weights so that cosine
        #  -weighed integral of f(x) over a hemisphere is
        #    F = SUM_i[ w_i * SUM_j( f(theta_i,phi_j) ) ]
        self.xquad_t, self.gq_w = gauleg_numpy(0, np.pi/2, self.nquad_t)
        self.wght_q = 2*self.gq_w*np.cos(self.xquad_t)*np.sin(self.xquad_t) / self.nquad_p

        self.rDDground = np.array( self.frtconf['rDDground'], dtype=float ) # ground reflectance, Dif->Dif (f77 rddgrou)
        self.rSDground = np.array( self.frtconf['rSDground'], dtype=float ) # ground reflectance, Sun->Dif (f77 rsdgrou)
        self.rDSground = np.array( self.frtconf['rDSground'], dtype=float ) # ground reflectance, Dif->Sensor (f77 rdsgrou)
        self.rSSground = np.array( self.frtconf['rSSground'], dtype=float  ) # ground reflectance, Sun->Sensor -- was "so" (Sun->Observer in f77 original code) (f77 rsogrou)
        self.SQratio = np.array( self.frtconf['SQratio'], dtype=float ) # f77 sqratio
        self.WaxRefrInd = np.array( self.frtconf['WaxRefrInd'], dtype=float )

        # mean LEAF, BRANCH, TRUNK REFLECTANCE AND TRANSMITTANCE
        self.EffRefrInd = [] # f77: rnlf
        self.LeafRefl = [] # not used in f77
        self.LeafReflLambert = [] # f77: rlfcl
        self.LeafTrans = [] # f77: tlfcl
        self.BranchRefl = [] # f77: rbrnc
        self.TrunkRefl = [] # f77: rtrnk
        for i,tc_dict in enumerate(self.frtconf["TreeClasses"]):
            wf = tc_dict['WaxCorrectionFactor'] if 'WaxCorrectionFactor' in tc_dict.keys() else 1
            self.EffRefrInd.append( wf*self.WaxRefrInd )
            self.EffRefrInd[-1][ np.nonzero( self.EffRefrInd[-1] < 1) ] = 1
            # setting WaxCorrectionFactor to 0 can thus be used to ignore specular reflectance
            #    (if WaxCorrectionFactor==1, refractive indices are equal and specular reflectance is 0)
            self.LeafRefl.append( tc_dict['LeafRefl'] ) # only adaxial reflectance used now
            self.LeafTrans.append( tc_dict['LeafTrans'] ) # f77: tlfcl
            # separate diffuse and specular reflectance components
            self.LeafReflLambert.append( self.LeafRefl[-1] -
                ((1 - self.EffRefrInd[-1])/(1 + self.EffRefrInd[-1]))**2 )
            self.BranchRefl.append( tc_dict['BranchRefl'] )
            self.TrunkRefl.append( tc_dict['TrunkRefl'] )
        # convert all lists of spectra to numpy arrays
        self.EffRefrInd = np.array( self.EffRefrInd, dtype=float )
        self.LeafRefl = np.array( self.LeafRefl, dtype=float )
        self.LeafReflLambert = np.array( self.LeafReflLambert, dtype=float )
        self.LeafTrans = np.array( self.LeafTrans, dtype=float )
        self.BranchRefl = np.array( self.BranchRefl, dtype=float )
        self.TrunkRefl = np.array( self.TrunkRefl, dtype=float )

        for i in range(self.ncl):
            if 'Description' not in self.frtconf['TreeClasses'][i].keys():
                self.frtconf['TreeClasses'][i]['Description'] = 'treeclass_'+str(i+1)
            if 'Age' not in self.frtconf['TreeClasses'][i].keys():
                self.frtconf['TreeClasses'][i]['Age'] = 0
        self.configuration_applied = True

    def configure_frt( self ):
        """ fill in quadrtatures, compute LAI and mean values, etc.
        """
        if not self.configuration_applied:
            # try to load the configuration, don't complain
            self.apply_config()
            if not self.frtconf_isread:
                print("configure_frt(): Error has occurred, aborting.")
                return
        ErrorsFixed = False # flag for signaling
        self.FGI = [] # Fisher's Grouping Index, 'glmp' in f77
        for i in range(self.ncl):
            #    FGI: see Eq. (8) by Nilson (1999).
            # The original limits in f77 code for FGI were 0.001 and 6.
            # Stick to these, although the former is never reached (TreeClumping=7)
            #     XXX give a warning somewhere if limits are reached
            if self.TreeClumping[i] < 0.358352:
                self.FGI.append( 6.0 )
                self.TreeClumping[i] = 0.358352
                print("\nToo small TreeClumping for tree class {:d}, limited to {:5.3f}".format(i,self.TreeClumping[i]), end="" )
                ErrorsFixed = True
            elif self.TreeClumping[i] > 6.9146699:
                self.FGI.append( 0.001 )
                self.TreeClumping[i] = 6.9146699
                print("\nToo large TreeClumping for tree class {:d}, limited to {:5.3f}".format(i,self.TreeClumping[i]), end="" )
                ErrorsFixed = True
            else:
                # brentq imported from scipy.optimize as the most "generic" method suggested there
                self.FGI.append ( brentq( CI_minfun, 0.0005, 7, args=self.TreeClumping[i] ) )
            # ShootLength cannot be zero, use a minimum of 1 cm.
            if self.correctFGI:
                # correct for the unrealistic case of canopy cover > crown cover by adjusting Fisher's index FGI
                # according to Nilson and Kuusk (2004, Agricultural and Forest Meteorology 124, 157–169)
                # if no overlapping crowns appears (canopy cover = crown cover), FGI = 1-CrownCover.
                # This sets the practical lower limit on the FGI regularity: even lower FGI would indicate even more
                # regular distribution, but because a limit has been reached, this would not decrease canopy transmittance.
                CrownCover_i = np.pi*(self.StandDensity[i]*self.CrownRadius[i]**2)
                if self.FGI[i] < (1-CrownCover_i):
                    self.FGI[i] = 1-CrownCover_i
                    print("\nCorrected FGI of tree class {:d} to {:5.3f}".format(i,self.FGI[i]), end="" )
                    ErrorsFixed = True
                    # If we change Fisher's Grouping Index, we need to change also tree distribution parameter TDP
                    if self.FGI[i]==1:
                        self.TreeClumping[i] = 1 #
                    else:
                        self.TreeClumping[i] = -np.log(self.FGI[i]) / (1-self.FGI[i])
                    #
                #
            if self.ShootLength[i] <= 0:
                self.ShootLength[i] = 0.01
                print("\nToo small ShootLength for tree class {:d}, set to {:5.3f}".format(i,self.ShootLength[i]), end="" )
                ErrorsFixed = True
            #

        # CALCULATE MEAN PARAMETERS FOR THE STAND (averaged over tree classes)
        self.strmean()
        #  strmean() sets up self.ulg, self.uuu, self.OpticalLAI, self.StandTreeDensity, self.MeanTreeHeight, self.MeanCrownLengthEllipsoid,
        #    self.MeanLengthCylinder, self.MeanCrownRadius, self.MeanDBH, self.MeanLeafMass, self.MeanSLW, self.CrownCover,
        #    self.CanopyCover, self.TrunkAreaBelowCrown, self.StandEffectiveLAI, self.StandLAI, self.StandBAI

        # gaps for the directions in the cubature
        gaps_q = hetk8s_singledirection( self.xquad_t, self.l_elli,
            self.ShootLength, self.StandDensity, self.TreeHeight, self.CrownLength1, self.CrownLength2, self.CrownRadius,
            self.DBH, self.FGI, self.ulg )

        # calculate also canopy recollision probability and effective (optical) LAI
        #    using Stenberg (2007) and Miller (1967), respectively
        self.DIFN = sum( self.wght_q*gaps_q )*self.nquad_p # DIFfuse Non-interceptance: DIFN=1-i0
        #  multiplication by nquad_p simulates the different azimuth angles in the quadrature,
        #               each giving the same result
        self.optPAI = - sum( self.wght_q*np.log(gaps_q) )*self.nquad_p

        # note: p is calculated using PAI, not just LAI
        self.p_Pola = 1-(1-self.DIFN)/(self.StandLAI+self.StandBAI)
        #           Stenberg 2007, Eq.17 [Remote Sensing of Environment 109, 221–224]

        # calculate gaps_s (gaps in solar direction)
        self.gaps_s = hetk8s_singledirection( [self.thets], self.l_elli,
            self.ShootLength, self.StandDensity, self.TreeHeight, self.CrownLength1, self.CrownLength2, self.CrownRadius,
            self.DBH, self.FGI, self.ulg )
        self.gaps_s = self.gaps_s[0]

        # spectral variables: truncate to [i1,i2]
        self.EffRefrInd = self.EffRefrInd[:,self.i1:self.i2+1]
        self.LeafRefl = self.LeafRefl[:,self.i1:self.i2+1]
        self.LeafReflLambert = self.LeafReflLambert[:,self.i1:self.i2+1]
        self.LeafTrans = self.LeafTrans[:,self.i1:self.i2+1]
        self.BranchRefl = self.BranchRefl[:,self.i1:self.i2+1]
        self.TrunkRefl = self.TrunkRefl[:,self.i1:self.i2+1]
        self.wl = self.wl[self.i1:self.i2+1]
        self.rDDground = self.rDDground[self.i1:self.i2+1] # ground reflectance, Dif->Dif (f77 rddgrou)
        self.rSDground = self.rSDground[self.i1:self.i2+1] # ground reflectance, Sun->Dif (f77 rsdgrou)
        self.rDSground = self.rDSground[self.i1:self.i2+1] # ground reflectance, Dif->Sensor (f77 rdsgrou)
        self.rSSground = self.rSSground[self.i1:self.i2+1] # ground reflectance, Sun->Sensor -- was "so" (Sun->Observer in f77 original code) (f77 rsogrou)
        self.SQratio = self.SQratio[self.i1:self.i2+1] # f77 sqratio
        self.nwl = len( self.wl )

        # Compute the mean values of optical parameters for computing diffuse
        # reflectance with the  2-stream submodel
        self.optmean()

        if ErrorsFixed:
            print("") # add a newline to output
        self.corrfact( ignore_configured=True ) # set the ignore_configured flag as configuration is loaded but flag not set yet
        self.frt_configured = True


    def reflectance( self, compute_gaps=True  ):
        """ compute forest directional reflectance and transmittance
        returns R and T, which are also stored in self

        if (not compute_gaps), canopy geometry is assumed unchanged, after the last computation,
            i.e., only changes in optical properties have been done
        """
        if not self.frt_configured:
            self.configure_frt()
            if not self.frt_configured:
                print("reflectance(): Error has occurred, aborting.")
                return
        if compute_gaps:
            # compute bidirectional gap probabilities
            self.bdgfu, self.bdgfd, self.btr1uk, self.btr1dk = hetk8s_integrated( [self.thetv], [self.phiv], self.thets, self.nzs,
                self.l_elli, self.ShootLength,  self.StandDensity, self.TreeHeight, self.CrownLength1, self.CrownLength2,
                self.CrownRadius, self.DBH, self.FGI, self.ulg, self.uuu )

            psgvuM = hetk8s_bidirectional( [self.thetv], [self.phiv], self.thets, self.nzs,
                self.l_elli, self.ShootLength, self.StandDensity, self.TreeHeight, self.CrownLength1, self.CrownLength2,
                self.CrownRadius, self.DBH, self.FGI, self.ulg )
            self.psgvu = psgvuM[0,0]

            self.gaps_vvec = hetk8s_singledirection( [self.thetv],
                self.l_elli, self.ShootLength, self.StandDensity, self.TreeHeight, self.CrownLength1, self.CrownLength2,
                self.CrownRadius, self.DBH, self.FGI, self.ulg)
            self.gaps_v = self.gaps_vvec[0]

        # Calculate single scattering from tree crowns
        bi0uM, bi0dM =  hetk8o( self.thets, [self.thetv], [self.phiv],
            self.StandDensity, self.uuu, self.bdgfu, self.bdgfd, self.btr1uk, self.btr1dk,
            self.pkhair, self.EffRefrInd,
            self.TrunkRefl, self.LeafBranchReflLamb, self.LeafBranchTrans )
        # bi0u: single scattering reflectance component (tree crowns); 3D ndarray (1x1xnwl)
        # bi0d: single scattering transmittance component (tree crowns) 3D ndarray (1x1xnwl)
        self.bi0u = bi0uM[0,0,:]
        self.bi0d = bi0dM[0,0,:]

        self.R1_c = self.SQratio*self.bi0u # 1st order reflectance component (in view direction) from crowns
        self.R1_g = self.SQratio* self.rSSground*self.psgvu #  1st order ground reflectance component
        self.R1 = self.R1_c + self.R1_g # total first-order reflectance component in view direction
        self.T1 = self.SQratio*self.bi0d #  1st order transmittance component

        self.rhd_hi_c, self.rhd_hi_g, self.thd_hi = twostr( self.thets, self.thetv, self.SQratio,
            self.OpticalLAI, self.TrunkAreaBelowCrown, self.CorrectionFactor,
            self.MeanRefl, self.MeanTrans, self.rDSground, self.rSDground, self.rDDground,
            self.gaps_s, self.gaps_v )

        # rhd_hi_c: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-directional),
        #     contribution by canopy, includes all-order diffuse-sky reflectance component
        # rhd_hi_g: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-directional),
        # contribution by ground (forest floor), includes all-order diffuse-sky reflectance component
        # thd_hi: higher-order component (excl. 1st-order direct) of canopy+soil hemispherical-directional transmittance
        #    includes all-order diffuse-sky reflectance component
        #    the higher-order components (*_hi) include all orders of reflectance of diffuse-sky radiation
        #  NOTE:  the higher-order components (*_hi) include
        #   * first-order reflectance of diffuse-sky radiation weighted by diffuse sky flux contribution
        #   *  higher order reflectance of direct sun and diffuse sky fluxes,
        #      weighted by direct and diffuse sky TOC flux contributions, respectively.

        self.R = self.R1 + self.rhd_hi_c + self.rhd_hi_g # total forest reflectance in view direction
        self.T = self.T1 + self.thd_hi  # total forest scattered transmittance in view direction

        return self.R, self.T


    def flux_reflectance( self, compute_gaps=True ):
        """ compute flux reflectance and transmittance
        returns R and T, which are also stored in self

        if (not compute_gaps), canopy geometry is assumed unchanged, after the last computation,
            i.e., only changes in optical properties have been done
        """
        if not self.frt_configured:
            self.configure_frt()
            if not self.frt_configured:
                print("flux_reflectance(): Error has occurred, aborting.")
                return
        if compute_gaps:
            self.precompute_quadrature_gaps()

        # Calculate single scattering from tree crowns
        bi0uQ, bi0dQ =  hetk8o( self.thets, self.xquad_t, self.xquad_p,
            self.StandDensity, self.uuu, self.bdgfuQ, self.bdgfdQ, self.btr1ukQ, self.btr1dkQ, self.pkhair, self.EffRefrInd,
            self.TrunkRefl, self.LeafBranchReflLamb, self.LeafBranchTrans )
        # bi0u: single scattering reflectance component (tree crowns); 3D ndarray (ntheta, nphi, nwl)
        # bi0d: single scattering transmittance component (tree crowns) 3D ndarray (ntheta, nphi, nwl)

        self.R1_cQ = self.SQratio*bi0uQ # 1st order reflectance component (in quadrature directions) from crowns
        self.R1_gQ = (self.SQratio* self.rSSground)[np.newaxis, np.newaxis, ...] * self.psgvuQ[..., np.newaxis] #  1st order ground reflectance component
        self.R1Q = self.R1_cQ + self.R1_gQ # total first-order reflectance component in quadrature directions
        self.T1Q = self.SQratio*bi0dQ #  1st order transmittance component

        # the diffuse fluxes (w/o 1st order component) are independent from azimuth
        #   so loop only over zenith angle
        self.rhd_hi_cQ = np.zeros( (self.nquad_t, self.nquad_p, self.nwl))
        self.rhd_hi_gQ = np.zeros( (self.nquad_t, self.nquad_p, self.nwl))
        self.thd_hiQ = np.zeros( (self.nquad_t, self.nquad_p, self.nwl))
        for i,(theta,gaps) in enumerate( zip(self.xquad_t,self.gaps_vvec) ):
            rhd_hi_c, rhd_hi_g, thd_hi = twostr(self.thets, theta, self.SQratio,
                self.OpticalLAI, self.TrunkAreaBelowCrown, self.CorrectionFactor,
                self.MeanRefl, self.MeanTrans, self.rDSground, self.rSDground, self.rDDground,
                self.gaps_s, gaps)
            self.rhd_hi_cQ[i,:,:] = rhd_hi_c
            self.rhd_hi_gQ[i,:,:] = rhd_hi_g
            self.thd_hiQ[i,:,:] =  thd_hi

        # rhd_hi_c: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-directional),
        # contribution by canopy, includes all-order diffuse-sky reflectance component
        # rhd_hi_g: higher-order reflectance (order higher than 1) in direction thetv (=hemispherical-directional),
        # contribution by ground (forest floor), includes all-order diffuse-sky reflectance component
        # thd_hi: higher-order component (excl. 1st-order direct) of canopy+soil hemispherical-directional transmittance
        #    includes all-order diffuse-sky reflectance component
        #    the higher-order components (*_hi) include all orders of reflectance of diffuse-sky radiation
        #  NOTE:  the higher-order components (*_hi) include
        #   * first-order reflectance of diffuse-sky radiation weighted by diffuse sky flux contribution
        #   *  higher order reflectance of direct sun and diffuse sky fluxes,
        #      weighted by direct and diffuse sky TOC flux contributions, respectively.

        self.RQ = self.R1Q + self.rhd_hi_cQ + self.rhd_hi_gQ # total forest reflectance in quadrature directions
        self.TQ = self.T1Q + self.thd_hiQ

        # integrate reflectance and transmittance
        self.R1_cF = ( self.wght_q[..., np.newaxis, np.newaxis]*self.R1_cQ ).sum(axis=0).sum(axis=0) # 1st order flux reflectance component from crowns
        self.R1_gF = ( self.wght_q[..., np.newaxis, np.newaxis]*self.R1_gQ ).sum(axis=0).sum(axis=0)  #  1st order flux reflectance component from ground
        self.Rhi_cF = ( self.wght_q[..., np.newaxis, np.newaxis]*self.rhd_hi_cQ ).sum(axis=0).sum(axis=0) # higher-order flux reflectance component from crowns
        self.Rhi_gF = ( self.wght_q[..., np.newaxis, np.newaxis]*self.rhd_hi_gQ ).sum(axis=0).sum(axis=0) # higher-order flux reflectance component from ground
        self.R1F = ( self.wght_q[..., np.newaxis, np.newaxis]*self.R1Q ).sum(axis=0).sum(axis=0) # total first-order flux reflectance component
        self.T1F = ( self.wght_q[..., np.newaxis, np.newaxis]*self.T1Q ).sum(axis=0).sum(axis=0)
        self.RF = ( self.wght_q[..., np.newaxis, np.newaxis]*self.RQ ).sum(axis=0).sum(axis=0)
        self.TF = ( self.wght_q[..., np.newaxis, np.newaxis]*self.TQ ).sum(axis=0).sum(axis=0)

        return self.RF, self.TF


    def precompute_quadrature_gaps( self ):
        """ precompute bidirectional gap probabilities for view and sun angles and quadrature
        this avoids calling hetk8s_* separately for each direction
        This is not required for simple bidirectional reflectance computations

        Note: called automatically when running self.flux_reflectance()
        """
        if not ( self.frt_configured ):
            self.configure_frt()

        self.gaps_vvec = hetk8s_singledirection( self.xquad_t,
            self.l_elli, self.ShootLength, self.StandDensity, self.TreeHeight, self.CrownLength1, self.CrownLength2,
            self.CrownRadius, self.DBH, self.FGI, self.ulg)
        self.psgvuQ = hetk8s_bidirectional( self.xquad_t, self.xquad_p, self.thets, self.nzs,
            self.l_elli, self.ShootLength, self.StandDensity, self.TreeHeight, self.CrownLength1, self.CrownLength2, self.CrownRadius, self.DBH, self.FGI, self.ulg )
        self.bdgfuQ, self.bdgfdQ, self.btr1ukQ, self.btr1dkQ = hetk8s_integrated( self.xquad_t, self.xquad_p, self.thets, self.nzs,
            self.l_elli, self.ShootLength, self.StandDensity, self.TreeHeight, self.CrownLength1, self.CrownLength2, self.CrownRadius, self.DBH, self.FGI, self.ulg, self.uuu )
        # self.bdgfuQF, self.bdgfdQF, self.btr1ukQF, self.btr1dkQF = hetk8s_integrated_f77( self.xquad_t, self.xquad_p, self.thets, self.nzs,
        #    self.l_elli, self.ShootLength, self.StandDensity, self.TreeHeight, self.CrownLength1, self.CrownLength2, self.CrownRadius, self.DBH, self.FGI, self.ulg, self.uuu )

    def corrfact( self, ignore_configured=False ):
        """ compute the spectral correction factor to conserve energy between the geometric-
            optic and two-stream submodels, self.CorrectionFactor

        See documentation and Mõttus, M., Stenberg, P. & Rautiainen, M. (2007). Photon recollision
        probability in heterogeneous forest canopies: compatibility with a hybrid GO model.
        Journal of Geophysical Research – Atmospheres 112, D03104, doi: 10.1029/2006JD007445.

        Fills the self.CorrectinFactor array

        Args:
        ignore_configured: whether to ignore the self.frt_configured flag. Used when calling from within self.configure_frt()
        """
        if not ignore_configured:
            if not self.frt_configured:
                self.configure_frt()
                # self.configure() calls corrfact, so do not recompute self.CorrectionFactor upon returning
                return

        if (self.wlcorr < -len(self.wl)) or (self.wlcorr > max(self.wl)):
            # no correction
            self.CorrectionFactor = np.ones_like( self.wl )
            return

        # ss is a very small number added so that leaf abedo would not be unity
        #   which causes instability in the twostream model implemented in layer.f
        ss = 1.0e-9
        # modify the needed values to get an energy-conservative forest
        # do this in a copy of self. Because all the modified values are arrays, a normal copy would suffice
        #   but it's safer to use deepcopy.
        F2 = deepcopy( self )
        F2.EffRefrInd = np.ones_like( self.EffRefrInd ) # eliminate specular reflectance
        # scale leaf reflectance and transmittance so that albedo -> 1. albedo=0 for which scaling doe not work
        w = self.LeafRefl + self.LeafTrans # leaf albedo
        i = np.nonzero( w == 0 )
        w[i] = 1
        reflfrac = self.LeafRefl / w # fraction of weighted reflectance to weighted total scattering, F = refl/(refl+tran)
        # set leaf R=reflfrac and T=1-reflfrac, so leaf W=1
        # make sure that leaves are not completely non-absorbing as this causes a singularity in layer.f
        F2.LeafRefl = reflfrac*(1-ss)
        F2.LeafTrans = (1.0-reflfrac)*(1-ss)
        F2.LeafRefl[i] = 0.5*(1-ss)
        F2.LeafTrans[i] = 0.5*(1-ss)
        F2.LeafReflLambert = F2.LeafRefl
        # WWW in the future, also leaf abaxial reflectance needs to be set
        # set branch and trunk R to 1
        F2.BranchRefl = np.ones_like(F2.LeafRefl)*(1-ss)
        F2.TrunkRefl = np.ones_like(F2.LeafRefl)*(1-ss)
        # set black ground and direct incidence; set cf to unity (no correction)
        F2.rDDground = np.zeros_like( F2.rDDground )
        F2.rSDground = np.zeros_like( F2.rSDground )
        F2.rDSground = np.zeros_like( F2.rDSground )
        F2.rSSground = np.zeros_like( F2.rSSground )
        F2.SQratio = np.ones_like( F2.SQratio )
        # set the i1 and i2 based on wlcorr
        if self.wlcorr < 0:
            F2.i1 = round( self.wlcorr ) - 1 # switch to python numbering
            F2.i2 = F2.i1 # only one wavelength used
        elif self.wlcorr > 0:
            # find closest wavelength
            F2.i1 = np.argmin( abs( self.wl-self.wlcorr) )
            F2.i2 = F2.i1 # only one wavelength used
        F2.wlcorr = -len(F2.wl)-100 # no correction for F2

        # configure F2: compute mean values etc.
        F2.configure_frt()
        F2.flux_reflectance()
        # calculate the correction factor
        ho_rt = F2.Rhi_cF + ( F2.TF - F2.T1F )
        if F2.i2 == F2.i1:
            # correction factor was computed for a single wavelength, but should be applied to all
            CF = ( 1 - (F2.R1_cF+F2.T1F+F2.gaps_s) ) / ho_rt
            self.CorrectionFactor = np.ones_like( self.wl )*CF[0]
        else:
            self.CorrectionFactor = ( 1 - (F2.R1_cF+F2.T1F+F2.gaps_s) ) / ho_rt


    def getarray( self, paramname ):
        """ Return the treeclass-specific variable paramname as an array.

        Missing values replaced by zero

        Args: paramname: string
        Returns: array of length self.ncl
        """
        return np.array( [t[paramname] if paramname in t.keys() else 0
            for t in self.frtconf["TreeClasses"]] )


    def strmean( self ):
        """ Calculate the mean values of structure parameters across tree classes

        based on the strmean subroutine in f77, A. KUUSK  23.09.1995

        fills in the following variables in self:
        ulg: uuu / 2 ??? ul * G  (f77: uuu)
        uuu: effective plant area density (leaves corrected for clumping + branches) (f77: uuu)
        OpticalLAI: LAI of a random canopy that causes extinction similar to that of the stand at thseff (= 40°) (f77: efflai)
        StandTreeDensity: total number of trees per m^2  (f77: sntr)
        MeanTreeHeight: mean tree height (f77: hmtree)
        MeanCrownLengthEllipsoid: mean length of the ellipsoid or conical part of crown (f77: vhekm)
        MeanLengthCylinder: mean length of the cylinder part of crown (f77: vhcilm)
        MeanCrownRadius: mean crown radius  (f77: rmcrown)
        MeanDBH: mean dbh  (f77: dbhmean)
        MeanLeafMass: mean leaf mass per m^2  (f77: rmassm)
        MeanSLW: mean SLW  (f77: slwm)
        CrownCover: crown cover (f77: vliit) a.a. crown cover.
           NOTE: CrownCover is defined as sum of crown areas over unit area, can be larger than one
        CanopyCover: canopy cover (f77: cano), ignoring within-crown gaps -- assuming opaque crowns
        MeanCrownVolume: mean crown volume NOTE: not used anywhere? (f77: ruum)
        TrunkAreaBelowCrown: total projected below-crown trunk area per m^2 calculated as pi*dbh*l0/2) (f77:tlty)
        StandEffectiveLAI: total effective LAI (corrected for shoot-level clumping) (f77: tlaief)
        StandLAI: total LAI (f77: tlai)
        StandBAI: total BAI (f77: tbai)
        StandEffectivePAI: total effective PAI, leaves + branches + stems  (f77: utot)
        """
        #                    (********** keskmiste arvutamine  **************)
        #                              mean structure parameters
        self.StandTreeDensity = 0
        self.MeanTreeHeight = 0
        self.MeanCrownLengthEllipsoid = 0
        self.MeanLengthCylinder = 0
        self.MeanCrownRadius = 0
        self.MeanDBH = 0
        self.MeanLeafMass = 0
        self.MeanSLW = 0
        self.CrownCover = 0
        self.CanopyCover = 0
        self.MeanCrownVolume = 0
        self.StandLAI = 0
        self.StandBAI = 0
        self.TrunkAreaBelowCrown = 0
        self.uuu = []
        self.ulg = []
        for i in range(self.ncl):
            stdi = self.StandDensity[i]
            self.StandTreeDensity +=  stdi # total number of trees per m^2
            ri2 = self.CrownRadius[i]**2
            vhi = self.CrownLength1[i] + self.CrownLength2[i] # total crown length
            self.MeanTreeHeight  += stdi*self.TreeHeight[i] # weighted sum of tree height
            self.MeanCrownLengthEllipsoid   += stdi*self.CrownLength1[i] # weighted sum of crown length (elli)
            self.MeanLengthCylinder  += stdi*self.CrownLength2[i] # weighted sum of crown length (cyl)
            if self.l_elli[i]:
                voi = 2*np.pi/3*ri2*self.CrownLength1[i] # crown volume
            else:
                voi = np.pi/30*ri2*(self.CrownLength1[i] + 3*self.CrownLength2[i]) # crown volume

            self.MeanCrownVolume += voi*stdi # total crown volume per m^2
            self.MeanCrownRadius += stdi*self.CrownRadius[i] # weighted sum of crown radii
            self.MeanDBH += stdi*self.DBH[i] # weighted sum of stem diam. [m]
            self.MeanLeafMass += stdi*self.DryLeafWeight[i] # leaf mass per m^2
            self.MeanSLW += stdi*self.SLW[i] # weighted sum of SLW
            r2stdi = stdi*ri2
            self.CrownCover += r2stdi # weighted sum of crown radii^2
            self.CanopyCover += self.TreeClumping[i]*r2stdi # used for calculating canopy closure
            self.TrunkAreaBelowCrown += stdi*self.DBH[i]*(self.TreeHeight[i] - vhi)
                                    # total projected trunk area below crown
            # the projection factor equals 1/2...
            rlaief = self.LAI[i]*self.ShootClumping[i] # effective LAI
            # self.StandEffectiveLAI += rlaief # total effective LAI
            self.StandLAI += self.LAI[i] # total LAI
            self.StandBAI += self.BAI[i] # total BAI
            claief = rlaief + self.BAI[i] # effective PAI
            self.uuu.append( claief/voi/stdi ) # effective leaf+branch area density
            self.ulg.append( self.uuu[i]/2 ) # uuu * G ?

        self.MeanTreeHeight /= self.StandTreeDensity
        self.MeanCrownLengthEllipsoid /= self.StandTreeDensity
        self.MeanLengthCylinder /= self.StandTreeDensity
        self.MeanCrownVolume /= self.StandTreeDensity
        self.MeanCrownRadius /= self.StandTreeDensity
        self.MeanDBH /= self.StandTreeDensity
        self.MeanLeafMass /= self.StandTreeDensity
        self.MeanSLW /= self.StandTreeDensity
        self.CrownCover *= np.pi
        self.CanopyCover = 1 - np.exp(-np.pi*self.CanopyCover)
        self.TrunkAreaBelowCrown *= np.pi*.5 #  projection factor 0.5

        self.StandEffectiveLAI  = sum( self.ShootClumping * self.LAI )
        self.StandEffectivePAI  = ( self.StandEffectiveLAI + sum(self.BAI)
            + self.TrunkAreaBelowCrown ) # total effective PAI, leaves + branches + stems

        # ********************  Calculation of LAI(eff) (09.2002) *****************
        thseff = 0.698 # view angle used = 40°
        gsf = 0.5 # gsf - G_spherical = 0.5
        cthets  = np.cos(thseff)
        tgths   = np.sin(thseff)/cthets
        alsc    = 0
        zz12    = 0

        for i in range(self.ncl):
            cellb  = self.CrownLength1[i]/2 # vertical semiaxis of crown ellipsoid
            vaheg  = 1 - self.FGI[i]
            clai = self.LAI[i]*self.ShootClumping[i] + self.BAI[i]  # effective LAI+BAI
            if self.l_elli[i]:
                szu1, vzui = pi11u(tgths, zz12, self.TreeHeight[i], cellb, self.CrownRadius[i])
            else:
                szu1, vzui = pi22u(tgths, zz12, self.TreeHeight[i], self.CrownLength1[i], self.CrownLength2[i], self.CrownRadius[i])
            if szu1 > 0:
            #                                    transmittance of a single crown
                a1gr = np.exp(-gsf*clai/(self.StandDensity[i]*szu1*cthets))
            else:
                a1gr = 1
            if np.abs(vaheg) > 0.1e-3:
            #                          parameter c in eq. (8) in Nilson (1999)
                b10i  = -np.log(1 - (1 - a1gr)*vaheg)/vaheg
            else:
                b10i  = 1 - a1gr
            alsc  += self.StandDensity[i]*szu1*b10i
        self.OpticalLAI =  cthets*alsc / gsf #  LAI of a random canopy that causes
        #                extinction similar to  that of the stand at thseff (= 40°)
        #  WWW  in the original version, the 2 (1/G) was missing
        return # =========================== strmean
    #

    def strmean_f77( self ):
        """ wrapper for the original fortran77 strmean. OBSOLETE. DO NOT USE
        Fills the same variables as strmean(), but appends 'F' to each name"""
        import strmean
        strmean.pidr.pi = np.pi
        strmean.pidr.dr = np.pi/180
        import enel3 # enel3 gives access to /frtpar/
        enel3.enel_incl()

        # f2py requires input arrays to be of exactly correct shape
        nclmax = enel3.frtpar.fnclmax[()]
        l_elli_F = np.empty(nclmax, dtype=bool)
        shl_F = np.F(nclmax)
        stdns_F = np.empty(nclmax)
        htr_F = np.empty(nclmax)
        hc1_F = np.empty(nclmax)
        hc2_F = np.empty(nclmax)
        rcr_F = np.empty(nclmax)
        dbh_F = np.empty(nclmax)
        rmass_F = np.empty(nclmax)
        slwcl_F = np.empty(nclmax)
        rlai_F = np.empty(nclmax)
        rbai_F = np.empty(nclmax)
        clmpst_F = np.empty(nclmax)
        clmpsh_F = np.empty(nclmax)
        glmp_F = np.empty(nclmax)

        ncl = len(self.l_elli)
        l_elli_F[0:ncl] = self.l_elli
        shl_F[0:ncl] = self.ShootLength
        stdns_F[0:ncl] = self.StandDensity
        htr_F[0:ncl] = self.TreeHeight
        hc1_F[0:ncl] = self.CrownLength1
        hc2_F[0:ncl] = self.CrownLength2
        rcr_F[0:ncl] = self.CrownRadius
        dbh_F[0:ncl] = self.DBH
        rmass_F[0:ncl] = self.DryLeafWeight
        slwcl_F[0:ncl] = self.SLW
        rlai_F[0:ncl] = self.LAI
        rbai_F[0:ncl] = self.BAI
        clmpst_F[0:ncl] = self.TreeClumping
        clmpsh_F[0:ncl] = self.ShootClumping
        glmp_F[0:ncl] = self.FGI

        self.ulgF, self.uuuF, self.OpticalLAIF, self.StandTreeDensityF, self.MeanTreeHeightF, self.MeanCrownLengthEllipsoidF, self.MeanLengthCylinderF, \
            self.MeanCrownRadiusF, self.MeanDBHF, self.MeanLeafMassF, self.MeanSLWF, self.CrownCoverF, self.CanopyCoverF, self.TrunkAreaBelowCrownF, \
            self.StandEffectiveLAIF, self.StandLAIF, self.StandBAIF = strmean.strmean( l_elli_F, ncl,
            stdns_F, htr_F, hc1_F, hc2_F, rcr_F, dbh_F, rmass_F, slwcl_F, rlai_F, rbai_F, clmpst_F, clmpsh_F, glmp_F)
        self.EffectivePAIF  = ( self.StandEffectiveLAIF + sum(self.BAI) + self.TrunkAreaBelowCrownF ) # total effective PAI, leaves + branches + stems
        # ===================== strmean_f77

    def optmean( self ):
        """the mean values of optical parameters    A. Kuusk 23.09.1995
        translation into Python: Matti Mõttus 2022

        defines the following variables:
        self.MeanLBRefl: average (area-weighted) reflectance of leaves+branches (i.e., excl. trunks) (f77: rleff)
        self.MeanLBTrans: average (area-weighted) tranmsmittance of leaves+branches (f77: tleff)
        self.MeanRefrIndex: average(area-weighted)  wax refractive index (f77: rneff)
        self.MeanTrunkRefl: average (area-weighted) trunk reflectance (f77: rty)
        self.LeafBranchReflLamb: area-weighted diffuse reflectace of a tree class, excl. trunks (for each tree class, f77: rrs)
        self.LeafBranchTrans: area-weighted diffuse transmittance of a tree class, excl. trunks (for each tree class, f77: ttt)
        self.MeanRefl: average (area-weighted) reflectance of canopy elements (f77: rteff)
        self.MeanTrans: average (area-weighted) leaf diffuse reflectance (ft: tteff) ?
        """

        if self.TrunkAreaBelowCrown != 0: # TrunkAreaBelowCrown: total projected below-crown trunk area per m^2
            sum2 = np.dot( 0.5*self.StandDensity*self.DBH*( self.TreeHeight-(self.CrownLength1+self.CrownLength2) ), self.TrunkRefl )
            self.MeanTrunkRefl = sum2*np.pi/self.TrunkAreaBelowCrown # average trunk reflectance
        else:
            self.MeanTrunkRefl = 0

        sumEffLAIBAI = self.StandEffectiveLAI + sum( self.BAI )
        EffLAI = self.ShootClumping*self.LAI
        self.LeafBranchReflLamb = []
        self.LeafBranchTrans = []
        for i in range(self.ncl):
            self.LeafBranchReflLamb.append( ( EffLAI[i]*self.LeafReflLambert[i] + self.BAI[i]*self.BranchRefl[i] ) \
                / (EffLAI[i]+self.BAI[i]) ) # area-weighted diffuse reflectace of a tree class, excl. trunks
            self.LeafBranchTrans.append( EffLAI[i]*self.LeafTrans[i]  / (EffLAI[i]+self.BAI[i]) )
        self.MeanLBRefl = ( np.dot(EffLAI,self.LeafRefl) + np.dot(self.BAI,self.BranchRefl) ) / sumEffLAIBAI # w/o trunks, old rleff
        self.MeanRefl = ( self.MeanLBRefl*sumEffLAIBAI + self.MeanTrunkRefl*self.TrunkAreaBelowCrown ) / self.StandEffectivePAI
        self.MeanTrans = np.dot(EffLAI,self.LeafTrans) / self.StandEffectivePAI # average diffuse leaf reflectance, f77: tteff
        self.MeanLBTrans = np.dot(EffLAI,self.LeafTrans) / sumEffLAIBAI # average diffuse transmittance of canopy elements
        self.LeafBranchReflLamb = np.array( self.LeafBranchReflLamb )
        self.LeafBranchTrans = np.array( self.LeafBranchTrans )

        if self.StandEffectiveLAI != 0:
            self.MeanRefrIndex = np.dot(EffLAI,self.EffRefrInd)/self.StandEffectiveLAI # average wax refractive index
        else:
            self.MeanRefrIndex = 1
        # in f77 version of rleff, specular component was added also to branches in MeanLBRefl. Recreate this for debugging
        self.rleff_old = np.dot( EffLAI+self.BAI, self.LeafBranchReflLamb) / sumEffLAIBAI + \
             ( (self.MeanRefrIndex - 1)/(self.MeanRefrIndex + 1))**2 # average reflectance of leaves+branches
        # self.rleff_old  should be equal to self.MeanReflLBF (available via self.optmean_f77())

    def optmean_f77( self ):
        """ wrapper for the original fortran77 optmean. OBSOLETE. DO NOT USE
        Fills the same variables as optmean(), but appends 'F' to each name
        """

        import optmean
        optmean.pidr.pi = np.pi
        optmean.pidr.dr = np.pi/180
        import enel3 # enel3 gives access to /frtpar/
        enel3.enel_incl()

        # f2py requires input arrays to be of exactly correct shape
        nclmax = enel3.frtpar.fnclmax[()]
        stdns_F = np.empty(nclmax)
        htr_F = np.empty(nclmax)
        hc1_F = np.empty(nclmax)
        hc2_F = np.empty(nclmax)
        rcr_F = np.empty(nclmax)
        dbh_F = np.empty(nclmax)
        rlai_F = np.empty(nclmax)
        rbai_F = np.empty(nclmax)
        clmpsh_F = np.empty(nclmax)

        ncl = len(self.l_elli)
        stdns_F[0:ncl] = self.StandDensity
        htr_F[0:ncl] = self.TreeHeight
        hc1_F[0:ncl] = self.CrownLength1
        hc2_F[0:ncl] = self.CrownLength2
        rcr_F[0:ncl] = self.CrownRadius
        dbh_F[0:ncl] = self.DBH
        rlai_F[0:ncl] = self.LAI
        rbai_F[0:ncl] = self.BAI
        clmpsh_F[0:ncl] = self.ShootClumping

        # output matrices
        Nwl = self.i2-self.i1+1
        self.MeanReflLBF = np.zeros( Nwl )
        self.MeanTransLBF = np.zeros( Nwl )
        self.MeanRefrIndexF = np.zeros( Nwl )
        self.MeanTrunkReflF= np.zeros( Nwl )
        self.ReflLambertLBF= np.zeros( (ncl, Nwl) )
        self.tttF = np.zeros( (ncl, Nwl) )
        self.MeanReflF = np.zeros( Nwl )
        self.MeanTransF = np.zeros( Nwl )

        rlfcl_i = np.empty(nclmax)
        tlfcl_i = np.empty(nclmax)
        rnlf_i = np.empty(nclmax)
        rbrnc_i = np.empty(nclmax)
        rtrnk_i = np.empty(nclmax)
        for i in range(self.i1,self.i2+1):
            rlfcl_i[0:ncl] = np.array(self.LeafReflLambert)[:,i]
            tlfcl_i[0:ncl] = np.array(self.LeafTrans)[:,i]
            rnlf_i[0:ncl] = np.array(self.EffRefrInd)[:,i]
            rbrnc_i[0:ncl] = np.array(self.BranchRefl)[:,i]
            rtrnk_i[0:ncl] = np.array(self.TrunkRefl)[:,i]

            rleff_i, tleff_i, rneff_i, rty_i, rrs_i, ttt_i, self.utot, rteff_i, tteff_i  = optmean.optmean(
            ncl, stdns_F, htr_F, hc1_F, hc2_F, dbh_F, rlai_F, rbai_F, clmpsh_F, self.TrunkAreaBelowCrown,
            self.StandEffectiveLAI, self.StandBAI, rlfcl_i, tlfcl_i, rnlf_i, rbrnc_i, rtrnk_i )

            self.MeanReflLBF[i] = rleff_i
            self.MeanTransLBF[i] = tleff_i
            self.MeanRefrIndexF[i] = rneff_i
            self.MeanTrunkReflF[i] = rty_i
            self.ReflLambertLBF[:,i] = rrs_i[0:ncl]
            self.tttF[:,i] = ttt_i[0:ncl]
            self.MeanReflF[i] = rteff_i
            self.MeanTransF[i] = tteff_i

