# FRT configuration is stored in a dictionary. This script will generate one in the dict frtconf
#   dicts can be stored as text files with pickle or json.dump() -- note: json.dump does not handle numpy arrays, which must be converted to lists first!
import numpy as np

frtconf = {} # create empty dict
# ---- technical parameters for the model
# Size of the cubature (quadrature) for flux computations. Optional, a "reasonable" cubature is predefined
frtconf['nquad_t'] = 7 # Optional: number of cubature knots oer the polar angle (theta)
frtconf['nquad_p'] = 7 # Optional: number of cubature knots oer the azimuth angle (phi) in the range [0,pi) (model is symmetrical)
frtconf['nzs'] = 4 # Optional: no. of canopy layers (used to integrate numerically 1st order reflectance). A reasonable default value has been set in code.
frtconf['wlcorr'] = 0 # Optional: the wavelength or thewavelength index for which the correction between first-order and higher order submodels is computed.
# By default, all wavelengths are used (value 0). A positive value indicates wavelength and is matched against wl (closest value is found); a negative value (should be integer) is index in wl.
# set to a very large negative value (out of bounds index, e.g., -1000000) to disable correction

frtconf['correctFGI'] = True # optional, defaults to True: whether to correct for the unrealistic case of canopy cover  > crown cover

# ---- environmand environmental parameters
frtconf['thetv'] = 0*np.pi/180 # view zenith angle [rad]
frtconf['thets'] = 40*np.pi/180 # solar zenith angle [rad]
frtconf['phiv'] = 0*np.pi/180 # relative view azimuth angle [rad]

# frtconf['wl'] = np.arange( 400, 2400, 5 )
# The wavelength used in the model. They serve as labels. The length of wl vector sets the
#   length requirement for all other input spectral vectors. For spectral data read from
#   external files, the values will be interpolated to the values in wl
frtconf['SQratio'] = np.ones_like( frtconf['wl'] )*1.0 # S/Q ratio: ratio of direct to total irradiance. Vector of the same length as wl

frtconf['wl'] = [ 450, 480, 510, 540, 570 ]
frtconf['i1'] = 1 # Optional: the the beginning of subset (first index) of the wavelengths used in the computations. Defaults to 1. Note: indexing starts from 1 (human, not python style)
frtconf['i2'] = len(frtconf['wl']) # optional: the the beginning of subset (first index) of the wavelengths used in the computations. Defaults to len(wl)


# ---- parameters quantifying the fores stand
frtconf['Description'] = "plot_813030"          # optional, not used in computations
frtconf['Age'] =  35                   # optional, not used in computations

# frtconf['WaxRefrInd'] = [1.4955, 1.491 , 1.4861, 1.4774, 1.4662]
# Wax refractive index (e.g., from the Prospect model) for the wavelngths in wl
frtconf['WaxRefrIndFile'] = "refrind.dat"
# If WaxRefrInd is not present in the dict, this file will be read. 2-column format, (wavelength, idx)

# Four separate spectra describe the forest floor. It may be reasonable to use the same spectrum for all (if no data available)
frtconf['rDDground'] = [0.0338, 0.0384, 0.0489, 0.0716, 0.0754] # Diffuse-Diffuse reflectance (=albedo). If not given, rDDgroundFile must be specified
frtconf['rDDgroundFile'] = "understory_fertility4.txt" # not used if rDDground is available
frtconf['rSDground'] = frtconf['rDDground'] # Sun-Diffuse reflectance (DHRF). If not given, rSDgroundFile must be specified
frtconf['rSDgroundFile'] = "understory_fertility4.txt" # not used if rSDground is available
frtconf['rDSground'] = frtconf['rDDground'] # Diffuse-Sensor reflectance (HDRF). If not given, rDSgroundFile must be specified
frtconf['rDSgroundFile'] = "understory_fertility4.txt" # not used if rDSground is available
frtconf['rSSground'] = frtconf['rDDground'] # Sun-Sensor reflectance (BDRF). If not given, rDDgroundFile must be specified
frtconf['rSSgroundFile'] = "understory_fertility4.txt" # not used if rSSground is available

frtconf['pkhair'] = 1 # Optional: leaf hair optical index (vaguely defined, set to unity by default)


# ---- individual tree class descriptions
frtconf['TreeClasses'] = [ {} ] # tree classes are dicts inside the list 'TreeClasses'
# The first tree class, python numbering starts from 0
frtconf['TreeClasses'][0]['Description'] = "pine" # arbitrary name
frtconf['TreeClasses'][0]['l_elli'] = True # Flag: wheter the class is ellipsoid (alternatively, cone+cylinder)
frtconf['TreeClasses'][0]['StandDensity'] =  0.0354 # stand density, m^-2
frtconf['TreeClasses'][0]['TreeHeight'] =  15.6 # tree height, m
frtconf['TreeClasses'][0]['CrownLength1'] = 8.8  # crown length, m -- ellipsoid or cone
frtconf['TreeClasses'][0]['CrownLength2'] = 0 # optional, crown length, cylinder (used only if l_elli == False)
frtconf['TreeClasses'][0]['CrownRadius'] = 2.12 # crown radius, m
frtconf['TreeClasses'][0]['DBH'] = 23.6  # trunk diameter (diameter at breast height), cm
frtconf['TreeClasses'][0]['DryLeafWeight'] = 9.3  # dry leaf weight per tree, kg
frtconf['TreeClasses'][0]['SLW'] = 161.3  # specific leaf weight, g m^-2. SLW and DryLeafWeight used to compute LAI
frtconf['TreeClasses'][0]['BAILAI'] = 0.18  # BAI / LAI [BAI = branch area index]
frtconf['TreeClasses'][0]['TreeClumping'] = 1.5  # tree distribution parameter, used to compute Fisher's Index. Small TreeClumping > large FI > clustered distribution
frtconf['TreeClasses'][0]['ShootClumping'] = 0.558  # shoot shading coefficient (ssc)
frtconf['TreeClasses'][0]['ShootLength'] =  0.1  # shoot length [m]. If no shoots, set to approx. leaf size
frtconf['TreeClasses'][0]['LeafRefl'] = [0.07636803, 0.07971808, 0.10241867, 0.17387372, 0.15803358] # if not given, LeafReflFile must be given
frtconf['TreeClasses'][0]['LeafTrans'] = [0.01553344, 0.02037519, 0.04619775, 0.12199694, 0.11064153] # if not given, LeafReflFile must be given
frtconf['TreeClasses'][0]['LeafReflFile'] = "leafspectrum_pine.txt" # not used if LeafRefl and LeafTrans are given.
#  Leaf spectral files are assumed to have wl in first column, followed by leaf reflectance and leaf transmittance
#     if a fourth column is present, the 2nd column is assumed leaf adaxial refl, and 4th column leaf abaxial reflectance
# frtconf['TreeClasses'][0]['BranchRefl'] = [0.3604, 0.3773, 0.3921, 0.4049, 0.4156] # if not given, BranchReflFile must be given
frtconf['TreeClasses'][0]['BranchReflFile'] = "branchspectrum_pine.txt" # branch reflectance file, optional: not used if BranchRefl given
# frtconf['TreeClasses'][0]['TrunkRefl'] = [0.3604, 0.3773, 0.3921, 0.4049, 0.4156] # if not given, TrunkReflFile must be given
frtconf['TreeClasses'][0]['TrunkReflFile'] = "trunkspectrum_pine.txt" # branch reflectance file, optional: not used if TrunkRefl given
frtconf['TreeClasses'][0]['WaxCorrectionFactor'] = 1  # Optional: wax refractive index correction factor. Default 1. Setting to 0 will remove specular component.
frtconf['TreeClasses'][0]['ScaleNeedle'] = True # whether to use ShootClumping and p-theory to scale needle albedo -> shoot albedo, relevant only if ssc<1 (conifers)

# second tree class, python #1
frtconf['TreeClasses'].append( {} )
frtconf['TreeClasses'][1]['Description'] = "spruce"
frtconf['TreeClasses'][1]['l_elli'] = False
frtconf['TreeClasses'][1]['StandDensity'] =  0.012
frtconf['TreeClasses'][1]['TreeHeight'] =  7.5
frtconf['TreeClasses'][1]['CrownLength1'] = 4.0
frtconf['TreeClasses'][1]['CrownLength2'] = 2.5
frtconf['TreeClasses'][1]['CrownRadius'] = 1.17
frtconf['TreeClasses'][1]['DBH'] = 8.3
frtconf['TreeClasses'][1]['DryLeafWeight'] = 2.38
frtconf['TreeClasses'][1]['SLW'] = 202.0
frtconf['TreeClasses'][1]['BAILAI'] = 0.18
frtconf['TreeClasses'][1]['TreeClumping'] = 1.5
frtconf['TreeClasses'][1]['ShootClumping'] = 0.644
frtconf['TreeClasses'][1]['ShootLength'] = 0.2
frtconf['TreeClasses'][1]['LeafReflFile'] = "leafspectrum_spruce.txt"
frtconf['TreeClasses'][1]['BranchReflFile'] = "branchspectrum_spruce.txt"
frtconf['TreeClasses'][1]['TrunkReflFile'] = "branchspectrum_spruce.txt"
frtconf['TreeClasses'][1]['ScaleNeedle'] = True

# third tree class, python #2
frtconf['TreeClasses'].append( {} )
frtconf['TreeClasses'][2]['Description'] = "birch"
frtconf['TreeClasses'][2]['l_elli'] = True
frtconf['TreeClasses'][2]['StandDensity'] =  0.0314
frtconf['TreeClasses'][2]['TreeHeight'] =  10.3
frtconf['TreeClasses'][2]['CrownLength1'] = 6.0
frtconf['TreeClasses'][2]['CrownRadius'] = 1.75
frtconf['TreeClasses'][2]['DBH'] = 11.4
frtconf['TreeClasses'][2]['DryLeafWeight'] = 1.21
frtconf['TreeClasses'][2]['SLW'] = 74.07
frtconf['TreeClasses'][2]['BAILAI'] = 0.15
frtconf['TreeClasses'][2]['TreeClumping'] = 1.5
frtconf['TreeClasses'][2]['ShootClumping'] = 1.0
frtconf['TreeClasses'][2]['ShootLength'] = 0.4
frtconf['TreeClasses'][2]['LeafReflFile'] = "leafspectrum_birch.txt"
frtconf['TreeClasses'][2]['BranchReflFile'] = "branchspectrum_birch.txt"
frtconf['TreeClasses'][2]['TrunkReflFile'] = "branchspectrum_birch.txt"
frtconf['TreeClasses'][2]['ScaleNeedle'] = False
