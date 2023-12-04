%.fpp : %.f90
	$(FPP) $< > $@

%.o : %.f90
	$(LD) -c $(F2FLAGS) $< -o $@

SRC_DIR  = ./
PY_SRC_DIR  = ../

QEDIR = $(or ${qedir}, ../../../)
include ${QEDIR}/make.inc
PYTHON := $(or ${PYTHON}, python3)

F2FLAGS = $(F90FLAGS) $(QEPY_INCS)

#export LDFLAGS := ${LDFLAGS}
# For the old version numpy (e.g 1.17.3)
export NPY_DISTUTILS_APPEND_FLAGS=1

ifeq ($(findstring -D__MPI, ${F2FLAGS}), -D__MPI)
   GOAL=mpi
else
   GOAL=serial
endif

ifeq ($(findstring ifort, ${F90}), ifort)
   FCOMPILER=intelem
else
   FCOMPILER=gnu95
endif

ifneq ($(findstring -std, ${CFLAGS}), -std)
   CFLAGS += -std=c99
endif

FPP = ${F90} -E -cpp $(F2FLAGS)
FSO = ${MPIF90} -shared $(CFLAGS) -Wl,-rpath,'$$ORIGIN/'
F2PY = CFLAGS='$(CFLAGS)' ${PYTHON} ../f90wrap2py.py \
				--fcompiler=$(FCOMPILER) --build-dir . \
				-c --f90exec=$(MPIF90) --f77exec=$(MPIF90) \
				--f90flags='${F2FLAGS}'

export LDFLAGS +=-Wl,-rpath,'$$ORIGIN/:$$ORIGIN/qepy'
$(info $(LDFLAGS))

WRAP_OBJS = $(QE_OBJS) $(QEPY_OBJS)

WRAP_FPP_FILES = $(WRAP_OBJS:%.o=%.fpp)

F90WRAP_FILES = f90wrap_*.f90

MODS_BASE = -lqepy_modules \
-lqepy_upflib \
-lqepy_xclib \
-lqepy_fftxlib \
-lqepy_laxlib \
-lqepy_utilxlib \
-lqepy_mbd

MODS_PW = -lqepy_pw -lqepy_ks_solvers -lqepy_dft_d3 $(MODS_BASE)

MODS_LR = -lqepy_lr_modules $(MODS_PW)

INCS_MODS_BASE = -I../modules \
-I../upflib \
-I../xclib \
-I../fftxlib \
-I../laxlib \
-I../utilxlib \
-I../mbd

INCS_MODS_PW = -I../pw -I../ks_solvers -I../dft_d3 $(INCS_MODS_BASE)

INCS_MODS_LR = -I../lr_modules $(INCS_MODS_PW)

INCS_MODS_BASE_QE = -I$(QEDIR)/upflib \
               -I$(QEDIR)/XClib \
               -I$(QEDIR)/Modules \
               -I$(QEDIR)/FFTXlib/src \
			   -I$(QEDIR)/LAXlib \
			   -I$(QEDIR)/UtilXlib \
			   -I$(QEDIR)/MBD \
			    $(FOX_MOD)

INCS_MODS_PW_QE = -I$(QEDIR)/PW/src -I$(QEDIR)/KS_Solvers -I$(QEDIR)/dft-d3 $(INCS_MODS_BASE_QE)

INCS_MODS_LR_QE = -I$(QEDIR)/LR_Modules $(INCS_MODS_PW_QE)

.PHONY: clean

default: all

clean:
	-rm -f ${F90WRAP_FILES} ${WRAP_FPP_FILES}
	-rm -f *.o *.mod
	-rm -rf src.* .libs .f2py_f2cmap lib${PNAME}.*so
	-rm -rf ../${PNAME} ../lib${PNAME}.*so 

all : ${F90WRAP_FILES}
	${F2PY} -m lib${PNAME} ${F90WRAP_FILES} $(QEPY_OBJS) $(LIB_OBJS) $(QEPY_LIBS) ${QELIBS}
	cp lib$(PNAME).*so ../lib$(PNAME).so
	if test -d $(PNAME); then mv $(PNAME) ../; fi
