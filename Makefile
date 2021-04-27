include ../make.inc

SRC_DIRS := ./src
PY_SRC_DIR := ./Python

QEMODS= ../PW/src/libpw.a ../Modules/libqemod.a ../KS_Solvers/libks_solvers.a ../FFTXlib/libqefft.a \
	   ../LAXlib/libqela.a ../UtilXlib/libutil.a ../dft-d3/libdftd3qe.a

QEINC =-I../Modules/ -I../KS_Solvers/ -I../FFTXlib/ \
	   -I../LAXlib/ -I../UtilXlib/ -I../dft-d3/ -I../PW/src/ -I.

MODULES_SOURCES = constants.f90 cell_base.f90 ions_base.f90
MODULES_FILES = $(addprefix ../Modules/,${MODULES_SOURCES})

PW_SOURCES = pwcom.f90 scf_mod.f90 read_file_new.f90 punch.f90
PW_FILES = $(addprefix ../PW/src/,${PW_SOURCES})

qepy_SOURCES= qepy_scatter_mod.f90 \
			  qepy_common.f90 qepy_mod.f90 \
			  qepy_setlocal.f90 qepy_v_of_rho.f90 qepy_pw2casino_write.f90 \
		      qepy_hinit1.f90 qepy_pwscf.f90 qepy_run_pwscf.f90 qepy_electrons.f90 \
			  qepy_forces.f90 qepy_stop_run.f90
qepy_FILES = $(addprefix ./src/,${qepy_SOURCES})
qepy_OBJS= $(qepy_SOURCES:%.f90=%.o)

WRAP_FILES = ${MODULES_FILES} ${PW_FILES} ${qepy_FILES}

F90WRAP_FILES = f90wrap_*.f90

WRAP_FPP_FILES = $(notdir $(WRAP_FILES:%.f90=%.fpp))

PWFLAGS = $(F90FLAGS) $(QEINC)
NOTI = -fPIC -nomodule -qopenmp -fpp -mcmodel=large
F2FLAGS = $(filter-out $(NOTI),$(PWFLAGS))

FPP = ${F90} -E -cpp $(F2FLAGS)

#$(info 'Install path :',${PY_INSTALL_DIR})

ifeq ($(PY_INSTALL_DIR), )
    PY2_DIR := $(shell python -m site --user-site)
    PY3_DIR := $(shell python3 -m site --user-site)
ifeq ($(PY3_DIR), )
    PY3_DIR = ${PY2_DIR}
else
    PY3_DIR = $(shell python3 -m site --user-site)
endif
else
    PY2_DIR = $(PY_INSTALL_DIR)
    PY3_DIR = $(PY_INSTALL_DIR)
endif

default: python

qepy_mod.o             : qepy_scatter_mod.o qepy_common.o
qepy_v_of_rho.o        : qepy_common.o
qepy_pw2casino_write.o : qepy_common.o
qepy_electrons.o       : qepy_common.o
qepy_pwscf.o           : qepy_common.o
qepy_hinit1.o          : qepy_setlocal.o


vpath %.f90 $(SRC_DIRS)

$(filter %.o,${qepy_OBJS}):%.o : %.f90
	$(LD) -c $(PWFLAGS) $< -o $@

${WRAP_FPP_FILES}: ${WRAP_FILES}
	for f in ${WRAP_FILES}; do $(FPP) $$f > $$(basename $${f%.*}).fpp; done

${F90WRAP_FILES}: ${qepy_OBJS} ${WRAP_FPP_FILES}
	f90wrap -v -m qepy ${WRAP_FPP_FILES} -k $(PY_SRC_DIR)/kind_map \
	    --init-file $(PY_SRC_DIR)/init.py -P

.PHONY: clean install python mpi python-clean python-install

uninstall: python-uninstall
install: python-install
clean: python-clean

python: ${F90WRAP_FILES}
	f2py-f90wrap --fcompiler=intelem --build-dir . \
		--opt=-O2 \
		-c -m _qepy ${F90WRAP_FILES} $(qepy_OBJS) \
		${F2FLAGS} ${QEMODS} $(LIBOBJS) $(QELIBS) -liomp5

mpi: ${F90WRAP_FILES}
	f2py-f90wrap \
		--fcompiler=intelem --build-dir . \
		--compiler=intelem \
		-c --f90exec=mpiifort --f77exec=mpiifort \
		--opt=-O2 \
		-m _qepy ${F90WRAP_FILES} $(qepy_OBJS) \
		${F2FLAGS} ${QEMODS} $(LIBOBJS) $(QELIBS)


python-install:
	cp -r qepy _qepy*.so ${PY2_DIR}
	cp -r qepy _qepy*.so ${PY3_DIR}

python-uninstall:
	-rm -rf ${PY2_DIR}/qepy
	-rm -rf ${PY3_DIR}/qepy
	-rm -rf ${PY2_DIR}/_qepy*.so
	-rm -rf ${PY3_DIR}/_qepy*.so

python-clean:
	-rm -f _qepy*.so ${F90WRAP_FILES} ${WRAP_FPP_FILES}
	-rm -rf qepy
	-rm -rf f90wrap_*.o qepy_*.o qepy_*.mod
	-rm -rf src.* .libs .f2py_f2cmap
