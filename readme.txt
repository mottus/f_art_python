readme.txt version: 2023/03/22

This folder contains the modified FRT (Forest Reflectance and Transmittance) code.
The following modifications have been made compared with the original version:

* Different submodels (6S, PROSPECT, LIBERTY, ACRMf) have been removed. Instead, FRT
  depends on externally provided (measured or modeled) spectra for forest floor,
  leaf reflectance and transmittance, and sky irradiance. 
* The input file format for the fortran version has been modified to allow 
  specification of the various spectra.
* A normalization has been applied to balance te first-order and higher-order
  scattering. See M�ttus, M., Stenberg, P. & Rautiainen, M. (2007), JGR � Atmospheres,
  doi:10.1029/2006JD007445 for details.
* The option to compute albedo, or flux reflectance, has been added: hence the name
  f-art.
* The model has been translated to python. The code includes fully functional 
  fortran code and various stages of integration of python and fortran via f2py
  - frtwrapper.py, allowing the call the main modules of FRT from python;
  - frtfunctions_f77.py, which, together with frtclass.py, allows running the program
      in python while using some of the fortran modules for computations;
  - frtcalss.py with frtfunctions_py.py is a full python translation of FRT. 
Computational efficiency decreases with the amount of python used in computations at 
the cost of flexibility, e.g., via using other models for generating leaf albedo, or 
substituting the two-stream submodel in FRT with a more flexible one. Alas, also the
probability of encountering a bug is larger in python.

To run anything but the pure python version, a fortran compiler (preferably gfortran)
is needed and the fortran version of FRT needs to be compiled with all the relevant 
object files (using the command 'make frt'). The relevant modules to be imported to 
python via f2py need also be compiled -- see the python code for the imported modules
and the fortran program files for compiling instructions. The hints there are given for
Windows. Under Windows, FRT has been compiled with gfortran installed using anaconda
(conda install -c conda-forge m2w64-gcc-fortran and conda install -c conda-forge make).
Various problems may arise when setting up an independent fortran compiler, e.g., 
modules involving input and output not loding in python etc. Compiling on Linux
should be easier.

Most of the documentation is provided in code. To start using the python code, open 
frtclass.py (the main code) and frtconf.py (configuraton data). The python version 
of FRT does not read any configuration files: it's configured using a dictionary
described in frtconf.py, which can be loaded and saved with the tools available in
python, such as json or pickle.

Some spectral files are given here for demonstration only. The origin of the data is
specifed in the files. Look up the original data for more spectra and information on
how they can be used.

Matti M�ttus (matti.mottus@gmail.com)

Excerpt from original FRT (pure fortran) readme.txt follows:
----------------------------------------------------------------------
Instructions and comments on the program frt

Andres Kuusk
Tartu Observatory
61602 T�ravere, Estonia
andres@aai.ee
www.aai.ee/~andres/
www.aai/bgf/


10 October 2002,
updated 29 July 2009,
updated 11 December 2012
______________________________________________________________________


To make the code extract source texts in a separate directory.

make frt
make clean removes object files,
make distclean removes object files and executables.

The code is debugged in the Linux environment (openSUSE-11.4,
GNU Fortran (SUSE Linux) 4.5.1 [gcc-4_5-branch revision 167585]. 
If you don't use the gfortran compiler then 
you have to modify the makefile.

The manual is available on-line:

http://www.aai.ee/
-> Structure and research projects
-> Remote sensing of vegetation
-> Projects
-> Forest reflectance model FRT 

Some details about the algorithms used may be found in the papers
Nilson, T. and Peterson, U., 1991. A forest canopy reflectance model
and a test case. Remote Sens. Environ. 37:131-142.
Kuusk, A. and Nilson, T., 2000. A directional multispectral forest
reflectance model. Remote Sens. Environ., 72:244-252.
Nilson, T. and Peterson, U., 1994. Age dependence of forest
reflectance: analysis of main driving factors. Remote Sens.
Environ., 48: 319-331.
