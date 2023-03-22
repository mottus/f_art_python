# FRT run from python: this program calls consecutively rd_cft, corrfact and comprt
# the latter is the module where frt computations take place
#  i.e., almost as good as calling frt from the command line, but with data available in python
# WARNING! USE AT YOUR OWN RISK! This is the least supported way to use FRT
# more preferred approaches
#  1) use the python-only version in frtclass.py
#  2) use the python version in frtclass.py, but modify it to import fortran modules via frtfunctions_f77
#  3) use the fortran-only FRT (via make frt)
# relies heavily on F2PY, https://www.numfys.net/howto/F2PY/
#  numpy can be called as, e.g.,
# f2py.exe -c --compiler=mingw32 -m rd_cfm rd_cfm.f
# I compiled the comprt module as
# f2py.exe -c --compiler=mingw32 -m comprt bck3.f bgrdd.o diffor.o enel3.o hetk8o.o hetk8s.o layer.o optmean.o pi11d.o pi11u.o pi22d.o pi22u.o rmsub.o spooi.o stem.o strmean.o twostr.o comprt.f
# more examples are at the top of individual .f files
# f2py can be called from python:
# numpy.f2py.compile(source, modulename='untitled', extra_args='', verbose=True, source_fn=None, extension='.f')
# Matti Mõttus 2022, 2023; based on fortran code by Andres Kuusk

import numpy as np
import os
import sys
import copy

configfilename = "infile_test.txt"
oufile = "outfile_test.txt"

# ftr_dir: the place where frt data and configuration files are
# frt_dir = '/homeappl/home/mottusma/frt/f-art'
frt_dir = "C:/data/koodid/f_art_test"
# frt_srcdir the place where the python and fortran modules are stored (and frt source code)
frt_srcdir = "C:/data/koodid/f_art_devel/f-art"

if not frt_dir in sys.path:
    sys.path.append(frt_dir)
if not frt_srcdir in sys.path:
    sys.path.append(frt_srcdir)

os.chdir(frt_dir) # frt needs to be run from its data folder

# import the fortran version of the frt model
#  three separate subroutines are called by the main function: rd_cfm, corrfact, and comprt
#    these are compiled and available as separate modules, although corrfact includes a full comprt
from frt_wrapper_functions import *
import rd_cfm
import corrfact
corrfact.pidr.pi = np.pi # pi and pidr are actually filled also in fortran code
corrfact.pidr.dr = np.pi/180
corrfact.volint.l_volint = 0 # to indicate that integration weights need to be computed
import comprt
comprt.volint.l_volint = 0 # to indicate that integration weights need to be computed

InputData=rd_cfm.xd_cfm( configfilename )
# subroutine xd_cfm( fname, XC, IC, SPCIN, DESC, lerr)
# InputData: ( XC, IC, SPCIN, DESC, lerr)
XC = InputData[0]
IC = InputData[1]
SPCIN = InputData[2]
DESC = InputData[3]
lerr = InputData[4]

DESC_2 = DESC.reshape(-1,77).view('S77')
DESC_3 = [ x[0].decode('utf8').strip() for x in DESC_2 ]

# example of reading /frtpar/ elements
nclmax = rd_cfm.frtpar.fnclmax[()] # zero-dimensional ndarray
nspchnl = rd_cfm.frtpar.fnspchnl[()]
# other variables fnclmax, fncub, fnspchnl, fnknotm, fnphim

ijob = IC[11+nclmax]
print( 'Input read, job ID='+str(ijob) )

#  =================================================================
#     Calculate correction factors

icorr = -1 # indicating no correction, set the actual value later
N_XOUT = 20
N_SPCOUT = 9
N_IC = 20+nclmax
XOUT = np.zeros( N_XOUT, order='F' )
SPCOUT = np.zeros( (nspchnl, N_SPCOUT), order='F' )
if (IC[1] == 0 ):
    # calculate correction for each wavelength separately: slow!
    print('Correction for each wavelength separately: slow!')
    icorr = IC[9] # the first wavelength
    CorrectionFactor = corrfact.corrfact( XC, IC, SPCIN, N_XOUT, N_SPCOUT )
    SPCIN[:,7+5*nclmax] = CorrectionFactor
    print('Correction factors calculated.')
elif ( IC[1] > IC[10] ):
    # the index is out of range
    print('No correction used.')
    # set correction factors to 1
    SPCIN[ IC[9]-1:IC[10], 7+5*nclmax ] = 1
    CorrectionFactor = np.ones( nspchnl )
else:
    # the wavelength index used for correction is stored in IC(2)
    print('Correction factors for wl #',IC[1] )
    IC_temp = copy.copy( IC )
    IC_temp[9] = IC[1]
    IC_temp[10] = IC[1]
    icorr = IC[1]
    CorrectionFactor = corrfact.corrfact( XC, IC_temp, SPCIN, N_XOUT, N_SPCOUT )
    print('Correction factors calculated.')
    SPCIN[:,7+5*nclmax] = CorrectionFactor[ IC[1]-1 ]
#  set flags to recompute everything
IC[3] = 1
IC[4] = 1
IC[5] = 1
IC[6] = 1
IC[7] = 1
IC[8] = 1

# create output file and store stand data, then call the comprt subroutine in fortran
print('Opening output file '+ oufile )
with open(oufile,'w') as f:
    print('%1s %11s %s' % ('#','name ', DESC_3[0] ), file=f )
    print('%1s %11s %i4' % ('#','age', round(XC[0,24]) ), file=f )
    if ( ijob < 2 ):
        # calling fortran functions, which have multidimensional variable-length inputs, does not work with f2py
        #   below is a workaround, which probably assumes that SPCOUT occupies continuous memory
        SPCOUT_reshape = np.reshape( SPCOUT, np.product(SPCOUT.shape) )
        comprt.comprt(XC, IC, SPCIN, XOUT, SPCOUT_reshape[0:nspchnl] )
        SPCOUT = np.reshape(SPCOUT_reshape,(N_SPCOUT,nspchnl)).transpose()
        # create some handy variables
        wl = SPCIN[ IC[9]-1:IC[10], 0]
        R = SPCOUT[ IC[9]-1:IC[10], 0]
        T = SPCOUT[ IC[9]-1:IC[10], 1]
        R1_c = SPCOUT[ IC[9]-1:IC[10], 2]
        R1_g = SPCOUT[ IC[9]-1:IC[10], 3]
        Rhi_c = SPCOUT[ IC[9]-1:IC[10], 4]
        Rhi_g = SPCOUT[ IC[9]-1:IC[10], 5]
        T1 = SPCOUT[ IC[9]-1:IC[10], 6]
        LeafReflEff = SPCOUT[ IC[9]-1:IC[10], 7]
        LeafTransEff =  SPCOUT[ IC[9]-1:IC[10], 8]
        # store stand structural data in output file
        print('%1s %11s %7.3f' % ('#','LAIeff_sh', XOUT[9]), file=f )
        print('%1s %11s %7.3f' % ('#','LAI', XOUT[10]), file=f )
        print('%1s %11s %7.3f' % ('#','PAI', XOUT[10]+XOUT[8]), file=f )
        print('%1s %11s %7.3f' % ('#','crown_cl', XOUT[11]), file=f )
        print('%1s %11s %7.3f' % ('#','can_cover', XOUT[12]), file=f )
        print('%1s %11s %7.3f' % ('#','DIFN', XOUT[13]), file=f )
        print('%1s %11s %7.3f' % ('#','PAI_opt', XOUT[14]), file=f )
        print('%1s %11s %7.3f' % ('#','p_Pola', XOUT[15]), file=f )
        print('%1s %11s %7.3f' % ('#','gap_Sun', XOUT[17]), file=f )
        print('%1s %11s %7.3f' % ('#','LAI40deg', XOUT[19]), file=f )
        print('%1s %11s %7.3f' % ('#','gap_view', XOUT[16]), file=f )
        print('%1s %11s %7.3f' % ('#','gap_bidi', XOUT[18]), file=f )
        if ( icorr != -1 ):
            print('%1s %7s %i %7.3f' % ('#','Corr@',int(SPCIN[icorr-1,0]), SPCIN[IC[9]-1,7+5*nclmax]), file=f )
        # write header line to output file
        print('%1s %4s %12s %12s %12s %12s %12s %12s %12s %8s %8s' %
            ('#','wl', 'HDRF', 'HDTF', 'R1_cr' , 'R1_gnd' , 'Rhi_cr' , 'Rhi_gnd', 'T1',
            'R_elem', 'T_elem'), file=f )
        # and, finally, store data
        for iwl in range(IC[9]-1,IC[10]):
            print('%6.1f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %8.3f %8.3f' %
                ( SPCIN[iwl,0], SPCOUT[iwl,0], SPCOUT[iwl,1], SPCOUT[iwl,2], SPCOUT[iwl,3],
                SPCOUT[iwl,4], SPCOUT[iwl,5], SPCOUT[iwl,6], SPCOUT[iwl,7], SPCOUT[iwl,8] ), file=f )
        #
    #
#