# Makefile to make frt.   A. Kuusk. 29.04.2005
#                         M. Mõttus 22.10.2020
# create system-agnostic RM https://stackoverflow.com/questions/4058840/makefile-that-distinguishes-between-windows-and-unix-like-systems
ifdef OS
   RM = del /Q
   FixPath = $(subst /,\,$1)
   F2PY =  f2py.exe -c --compiler=mingw32
else
   ifeq ($(shell uname), Linux)
      RM = rm -f
      FixPath = $1
	  F2PY =  f2py -c 
   endif
endif

GG = gfortran
# GG     = g77

# GFLAGS = -c -O0 -Wunused -Wuninitialized -I.
# oli enne:
# GFLAGS = -c -O2 -Wunused -Wuninitialized -I.
# oli enne: GFLAGS = -g -c  -I.
GFLAGS = -g -c -Wunused -Wuninitialized  -fPIC -I.
# -fPIC may be required for linking with python
LFLAGS = -o
# LFLAGS = -static -o
# GFLAGS = -c -g -O0 -Wunused -I.
# LFLAGS = -g -o

###############

OBJ1  = frt.o comprt.o
OBJ   = bck3.o bgrdd.o corrfact.o \
 enel3.o hetk8o.o hetk8s.o \
 layer.o optmean.o rmsub.o \
 spooi.o strmean.o \
 rd_cfm.o twostr.o 

INCL = frtpar.h
INCL2 = cf_new.h

frt.o: frt.f $(INCL)
	$(GG) $(GFLAGS) $<

%.o: %.f $(INCL)
	$(GG) $(GFLAGS) $<

frt: $(OBJ) $(OBJ1) $(INCL) $(INCL2)
	$(GG) $(LFLAGS) $@ $(OBJ1) $(OBJ)

#frt: $(OBJ) $(OBJ1) $(INCL)
#	$(GG) $(LFLAGS) $@ $(OBJ1) $(OBJ); \
#	chmod g+rx frt

xd_cfm.pyd: rd_cfm.o
	$(F2PY) -m xd_cfm rd_cfm.o xd_cfm.f

spooi.pyd: rmsub.o
	$(F2PY) -m spooi rmsub.o spooi.f
	
enel3.pyd: spooi.o bck3.o rmsub.o 
	$(F2PY) -m enel3 spooi.o bck3.o rmsub.o enel3.f
	
bck3.pyd: spooi.o
	$(F2PY)  -m bck3 rmsub.o spooi.o bck3.f
	
pythonlibs: xd_cfm.pyd spooi.pyd enel3.pyd bck3.pyd

clean:
	$(RM) *.o core *.pyd

distclean:
	$(RM) *.o frt frt.exe core *.pyd
