# Makefile to make frt.   A. Kuusk. 29.04.2005
#                         M. Mõttus 22.10.2020
# create system-agnostic RM https://stackoverflow.com/questions/4058840/makefile-that-distinguishes-between-windows-and-unix-like-systems
ifdef OS
   RM = del /Q
   FixPath = $(subst /,\,$1)
else
   ifeq ($(shell uname), Linux)
      RM = rm -f
      FixPath = $1
   endif
endif

SHELL  = /bin/sh

GG     = gfortran
# GG     = g77

# GFLAGS = -c -O0 -Wunused -Wuninitialized -I.
# oli enne:
# GFLAGS = -c -O2 -Wunused -Wuninitialized -I.
# oli enne: GFLAGS = -g -c  -I.
GFLAGS = -g -c -Wunused -Wuninitialized  -fPIC -I.
# -fPIC on vajalik pythoniga linkimiseks
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
 
INCL1 = frtpar.h

frt.o: frt.f $(INCL1)
	$(GG) $(GFLAGS) $<

%.o: %.f $(INCL) $(INCL1)
	$(GG) $(GFLAGS) $<

frt: $(OBJ) $(OBJ1) $(INCL) 
	$(GG) $(LFLAGS) $@ $(OBJ1) $(OBJ)

#frt: $(OBJ) $(OBJ1) $(INCL)
#	$(GG) $(LFLAGS) $@ $(OBJ1) $(OBJ); \
#	chmod g+rx frt

clean:
	$(RM) *.o core *.pyd

distclean:
	$(RM) *.o frt frt.exe core *.pyd
