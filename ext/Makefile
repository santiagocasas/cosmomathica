# Define your C compiler
CC = gcc -O2 -Wl,--no-as-needed -Wall -lgfortran -fPIC -fopenmp -lm -ldl 

# C++ compiler for copter
CXX = g++ -O2 -Wl,--no-as-needed -fopenmp -Wall -pedantic -lm -ldl

# Define your Fortran compiler (if you want to use CAMB) Must match CAMB
# compiler
F90C = gfortran -O3 -fopenmp -ffast-math -fmax-errors=4 -ffixed-line-length-0 -std=legacy

include Makefile_mathlink 
# No changes necessary below this line
# ====================================

CADDSDIR = ${MLINKDIR}/${SYS}/CompilerAdditions
INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}
MPREP = ${CADDSDIR}/mprep

OFILES = math_link.o math_link_tm.o 


# Check if the external programs are present

ifeq ($(wildcard tf/tf_fit.c),)
  TRANSFER = 
else
  TRANSFER = tf
  OFILES += $(TRANSFER)/tf_fit.o $(TRANSFER)/power.o
  DFLAGS += -DTRANSFER
endif 

ifeq ($(wildcard halofit/smith2.c),)
  HALOFIT = 
else
  HALOFIT = halofit
  OFILES += $(HALOFIT)/smith2.o 
  DFLAGS += -DHALOFIT
endif

ifeq ($(wildcard CosmicEmulator/emu.c),)
  COSMICEMU = 
else
  COSMICEMU = CosmicEmulator
  OFILES += $(COSMICEMU)/hubble.o $(COSMICEMU)/emu.o 
  EMULIBS = -lgsl -lgslcblas 
  DFLAGS += -DCOSMICEMU
endif

ifeq ($(wildcard FrankenEmu/emu.c),)
  FRANKENEMU = 
else
  FRANKENEMU = FrankenEmu
  FRANKENOFILES += math_link2.o math_link_tm.o $(addprefix $(FRANKENEMU)/, hubble.o emu.o emu_noh.o)
  EMULIBS = -lgsl -lgslcblas 
endif

ifeq ($(wildcard camb/camb.f90),)
  CAMB = 
else
  CAMB = camb
  OFILES += camb_wrapper.o 
  CAMBOFILES = $(CAMB)/Release/*.o
  DFLAGS += -DCAMB
endif

ifeq ($(wildcard class/Makefile),)
  CLASS = 
else
  CLASS = class
  DFLAGS += -DCLASS
  OFILES += class_wrapper.o  #$(CLASS)/build/output.o
  EXTRALIBS += $(CLASS)/libclass.a
endif

all: math_link math_link2

math_link: $(HALOFIT) $(COSMICEMU) $(CAMB) $(CLASS) $(OFILES) Makefile
	${CC} -DMLINTERFACE=4 ${DFLAGS} -I${INCDIR} -L${LIBDIR} -l${MLLIB} -luuid \
		$(OFILES) $(CAMBOFILES) -o $@ $(EXTRALIBS)  \
		$(EMULIBS) $(COPLIBS)

# The second math_link is for FrankenEmu only, since it contains symbols
# with the same name as symbols from CosmicEmu
math_link2: math_link2.o $(FRANKENEMU) $(FRANKENOFILES) Makefile
	${CC} -DMLINTERFACE=4 -I${INCDIR} -L${LIBDIR} -l${MLLIB} -luuid $(FRANKENOFILES) \
			-DFRANKENEMU -o $@ $(EXTRALIBS) $(EMULIBS)

math_link2.o: math_link.c
	${CC} -DFRANKENEMU -c -I$(INCDIR) $< -o $@

math_link_tm.c: math_link.tm
	${MPREP} $? -o $@


# TRANSFER is built here
%.o: %.c
	$(CC) $(DFLAGS)  -c -I$(INCDIR) $< -o $@

camb_wrapper.o: camb_wrapper.f90 Makefile
	$(F90C) -c $< -o $@ -I$(CAMB)/Release 

class_wrapper.o: class_wrapper.c Makefile
	$(CC) -c -I$(CLASS)/include $< -o $@

$(HALOFIT):
	+$(MAKE) -C $@

$(COSMICEMU):
	+$(MAKE) -C $@

$(CAMB):
	+$(MAKE) -C $@
#	this script changes the ifort compiler to the gfortran compiler. Works with 2013 versions of CAMB
#       ./patchmakefile.sh

$(CLASS):
	+$(MAKE) -C $@


download:
	./download.sh


clean:
	rm -f *mod *o *tm.c math_link
	if [ ! -z $(TRANSFER) ]; then rm -f $(TRANSFER)/*o; fi
	rm -f $(COSMICEMU)/*{o,exe}
	if [ ! -z $(FRANKENEMU) ]; then $(MAKE) -C $(FRANKENEMU) clean ; fi
	if [ ! -z $(HALOFIT) ]; then $(MAKE) -C $(HALOFIT) clean ; fi
	if [ ! -z $(CAMBROOT) ]; then $(MAKE) -C $(CAMBROOT) clean ; fi
	if [ ! -z $(CLASS) ]; then $(MAKE) -C $(CLASS) clean ; fi

.PHONY: clean $(HALOFIT) $(COSMICEMU) $(CAMB) $(CLASS) download 
