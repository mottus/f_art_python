import os
import sys
# the following import are required to get the location of this file
from inspect import getsourcefile
from os.path import abspath

pythonconfname = "frtconf.py" # python script to create the configuration dict for option A below
configfilename = "in_frt_demo" # fortran-style "new" config file for option B

current_dir = os.path.dirname( abspath(getsourcefile(lambda:0)) )
# frt_srcdir is the location of python modules and functions (both compiled and not)
#  it's needed to import the python and pre-compiled fortran functions
frt_srcdir = os.path.join( current_dir, "f_art_python" )
if frt_srcdir not in sys.path:
    sys.path.append(frt_srcdir)

from frtclass import frt_model

frt_datadir = os.path.join( current_dir, "data" ) # the place where the FRT data files are (spectra etc.)

# In all situations, a FRT model instance needs to be generated
G = frt_model() # if no fortran is used, no arguments to frt_srcdir are required

# OPTION A
# Load the frtconf dict containing model configuration
# Here,it's defined in a python script in the same folder as this file
# exec( open( os.path.join(current_dir,"frtconf.py") ).read() )
# G.load_conf( frtconf , frt_datadir )
# G.configure_frt() # this is not needed as configure_frt called automatically by reflectance()
#print("computing reflectance for G...", end="")
# G.reflectance() # compute reflectance of the forest described in F
# print(" done. see G.R and G.T")

# OPTION B: use Fortran77 FRT "new" input file (text file)
# this requires fortran and compiled frt rd_cfm.f, xd_cfm.f

H = frt_model( frt_srcdir ) # if fortran is used, frt_srcdir needs to be given
H.read_conf( os.path.join( current_dir, configfilename ), frt_datadir )
H.configure_frt()
H.reflectance()