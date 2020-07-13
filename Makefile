include ../../make.inc

SRC_DIRS := ./src
PY_SRC_DIR := ./Python

QEMODS= ../src/libpw.a ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a ../../FFTXlib/libqefft.a \
	   ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a

QEINC =-I../../Modules/ -I../../KS_Solvers/ -I../../FFTXlib/ \
	   -I../../LAXlib/ -I../../UtilXlib/ -I../../dft-d3/ -I../src/ -I.

MODULES_SOURCES = read_cards.f90 read_input.f90 read_namelists.f90 read_pseudo.f90
MODULES_FILES = $(addprefix ../../Modules/,${MODULES_SOURCES})

PW_SOURCES = pwcom.f90 scf_mod.f90 electrons.f90
PW_FILES = $(addprefix ../src/,${PW_SOURCES})

PWPY_SOURCES= pwpy_electrons.f90 pwpy_forces.f90 pwpy_mod.f90 pwpy_pw2casino_write.f90 pwpy_pwscf.f90 pwpy_run_pwscf.f90 pwpy_sum_band.f90 pwpy_v_of_rho.f90
PWPY_FILES = $(addprefix ./src/,${PWPY_SOURCES})

PWPY_OBJS= $(addsuffix .o, $(basename ${PWPY_SOURCES}))

WRAP_SOURCES = ${MODULES_SOURCES} ${PW_SOURCES} ${PWPY_SOURCES}
WRAP_FILES = ${MODULES_FILES} ${PW_FILES} ${PWPY_FILES}

F90WRAP_FILES = f90wrap_*.f90

PWFLAGS = $(F90FLAGS) $(QEINC)
NOTI = -fPIC -nomodule -qopenmp -fpp -mcmodel=large
F2FLAGS = $(filter-out $(NOTI),$(PWFLAGS))

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


vpath %.f90 $(SRC_DIRS)

default: python

$(filter %.o,${PWPY_OBJS}):%.o : %.f90
	$(LD) -c $(PWFLAGS) $< -o $@

${F90WRAP_FILES}: ${PWPY_OBJS}
	f90wrap -v -m pwscfpy ${WRAP_FILES} -k $(PY_SRC_DIR)/kind_map \
	    --init-file $(PY_SRC_DIR)/init.py

.PHONY: clean install python python-clean python-install

install: python-install
clean: python-clean

python: ${F90WRAP_FILES}
	f2py-f90wrap --fcompiler=intelem --build-dir . \
		-c -m _pwscfpy ${F90WRAP_FILES} $(PWPY_OBJS) \
		${F2FLAGS} ${QEMODS} $(LIBOBJS) $(QELIBS) -liomp5

python-install:
	cp -r pwscfpy.py _pwscfpy*.so ${PY2_DIR}
	cp -r pwscfpy.py _pwscfpy*.so ${PY3_DIR}

python-clean:
	-rm -f _pwscfpy*.so ${F90WRAP_FILES}
	-rm -rf pwscfpy pwscfpy.py
	-rm -rf f90wrap_*.o pwpy_*.o pwpy_*.mod
	-rm -rf src.* .libs .f2py_f2cmap

