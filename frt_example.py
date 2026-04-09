from pathlib import Path
from f_art.frtclass import frt_model


current_dir = Path(__file__).resolve().parent


pythonconfname = "frtconf.py" # python script to create the configuration dict for option A below
configfilename = "in_frt_demo" # fortran-style "new" config file for option B

# frt_srcdir is the location of python modules and functions (both compiled and not)
#  it's needed to import the python and pre-compiled fortran functions
frt_srcdir = str( current_dir / 'f_art' )

# frt_datadir: the place where the FRT data files are (spectra etc.)
frt_datadir = str( current_dir / 'data' )

# frt_sampleconfigdir holds sample json files with frt inputs for testing
frt_sampleconfigdir = str( current_dir / 'test_forests' )
# a folder containing additional (external) spectra
spectradir = str( Path(r'C:\Users\mmattim\2026\kasikirjad\manuscript_FisherIndex\f_art_python\data') )

# In all situations, a FRT model instance needs to be generated
G = frt_model() # if no fortran is used, no arguments to frt_srcdir are required

# OPTION A
# Load the frtconf dict containing model configuration
# Here,it's defined in a python script in the same folder as this file
exec( open( str(current_dir / "frtconf.py") ).read() )
G.load_conf( frtconf , frt_datadir )
G.configure_frt() # this is not needed as configure_frt called automatically by reflectance()
print("computing reflectance for G...", end="")
G.reflectance() # compute reflectance of the forest described in F
print(" done. see G.R and G.T")

print("\n --- using json file from "+frt_sampleconfigdir)
G2 = frt_model( frt_datadir=spectradir, frtconffile=frt_sampleconfigdir+'/a_test.json' )
G2.reflectance()
print(" See G2.R and G2.T")

# OPTION B: use Fortran77 FRT "new" input file (text file)
# this requires fortran and compiled frt rd_cfm.f, xd_cfm.f
# H = frt_model( frt_srcdir ) # if fortran is used, frt_srcdir needs to be given
# H.read_conf( os.path.join( current_dir, configfilename ), frt_datadir )
# H.configure_frt()
# H.reflectance()