from frtclass import frt_model
import os
import sys


configfilename = "infile_test.txt"

current_dir = os.path.dirname(__file__)
frt_dir = os.path.join( current_dir, "data" ) # the place where the FRT data files are (spectra etc.)

# reading a configuration file, requires fortran and compiled frt
# F =frt_model( frt_srcdir  )
# F.read_conf( os.path.join(frt_dir,configfilename) )
# print("computing reflectance for F...", end="")
# F.reflectance()
# print(" done. see F.R and F.T")


# First, a model needs to be generated
G = frt_model( None )

# frt_srcdir is needed if compiled fortran modules are used
# frt_srcdir = current_dir # the place where the compiled f77 modules are, assuming in the same folder as this script
# H = frt_model( frt_srcdir )
# H.read_conf( os.path.join( current_dir, "in_new_demo" ) )

# Load the frtconf dict containing model configuration
# We assume that it is in the same folder as this file
exec( open( os.path.join(current_dir,"frtconf.py") ).read() )

G.load_conf( frtconf , frt_dir )
# G.configure_frt() # this is not needed as configure_frt called automatically by reflectance()
print("computing reflectance for G...", end="")
G.reflectance() # compute reflectance of the forest described in F
print(" done. see G.R and G.T")