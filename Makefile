include ../../make.inc

SRC_DIRS := ./src
PY_SRC_DIR := ./Python

QEMODS= ../src/libpw.a ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a ../../FFTXlib/libqefft.a \
	   ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a

QEINC =-I../../Modules/ -I../../KS_Solvers/ -I../../FFTXlib/ \
	   -I../../LAXlib/ -I../../UtilXlib/ -I../../dft-d3/ -I../src/ -I.

MODULES_SOURCES = constants.f90 cell_base.f90 ions_base.f90
MODULES_FILES = $(addprefix ../../Modules/,${MODULES_SOURCES})

PW_SOURCES = pwcom.f90 scf_mod.f90 read_file_new.f90
PW_FILES = $(addprefix ../src/,${PW_SOURCES})

PWPY_SOURCES= pwpy_scatter_mod.f90 \
			  pwpy_common.f90 pwpy_mod.f90 \
			  pwpy_setlocal.f90 pwpy_v_of_rho.f90 pwpy_pw2casino_write.f90 \
		      pwpy_hinit1.f90 pwpy_pwscf.f90 pwpy_run_pwscf.f90 pwpy_electrons.f90 \
			  pwpy_forces.f90 pwpy_stop_run.f90
PWPY_FILES = $(addprefix ./src/,${PWPY_SOURCES})
PWPY_OBJS= $(PWPY_SOURCES:%.f90=%.o)

WRAP_FILES = ${MODULES_FILES} ${PW_FILES} ${PWPY_FILES}

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

pwpy_mod.o             : pwpy_scatter_mod.o pwpy_common.o
pwpy_v_of_rho.o        : pwpy_common.o
pwpy_pw2casino_write.o : pwpy_common.o
pwpy_electrons.o       : pwpy_common.o
pwpy_pwscf.o           : pwpy_common.o
pwpy_hinit1.o          : pwpy_setlocal.o


vpath %.f90 $(SRC_DIRS)

$(filter %.o,${PWPY_OBJS}):%.o : %.f90
	$(LD) -c $(PWFLAGS) $< -o $@

${WRAP_FPP_FILES}: ${WRAP_FILES}
	for f in ${WRAP_FILES}; do $(FPP) $$f > $$(basename $${f%.*}).fpp; done

${F90WRAP_FILES}: ${PWPY_OBJS} ${WRAP_FPP_FILES}
	f90wrap -v -m pwscfpy ${WRAP_FPP_FILES} -k $(PY_SRC_DIR)/kind_map \
	    --init-file $(PY_SRC_DIR)/init.py -P

.PHONY: clean install python mpi python-clean python-install

uninstall: python-uninstall
install: python-install
clean: python-clean

python: ${F90WRAP_FILES}
	f2py-f90wrap --fcompiler=intelem --build-dir . \
		--opt=-O2 \
		-c -m _pwscfpy ${F90WRAP_FILES} $(PWPY_OBJS) \
		${F2FLAGS} ${QEMODS} $(LIBOBJS) $(QELIBS) -liomp5

mpi: ${F90WRAP_FILES}
	f2py-f90wrap \
		--fcompiler=intelem --build-dir . \
		--compiler=intelem \
		-c --f90exec=mpiifort --f77exec=mpiifort \
		--opt=-O2 \
		-m _pwscfpy ${F90WRAP_FILES} $(PWPY_OBJS) \
		${F2FLAGS} ${QEMODS} $(LIBOBJS) $(QELIBS)


python-install:
	cp -r pwscfpy _pwscfpy*.so ${PY2_DIR}
	cp -r pwscfpy _pwscfpy*.so ${PY3_DIR}

python-uninstall:
	-rm -rf ${PY2_DIR}/pwscfpy
	-rm -rf ${PY3_DIR}/pwscfpy
	-rm -rf ${PY2_DIR}/_pwscfpy*.so
	-rm -rf ${PY3_DIR}/_pwscfpy*.so

python-clean:
	-rm -f _pwscfpy*.so ${F90WRAP_FILES} ${WRAP_FPP_FILES}
	-rm -rf pwscfpy
	-rm -rf f90wrap_*.o pwpy_*.o pwpy_*.mod
	-rm -rf src.* .libs .f2py_f2cmap
