include make.qe.inc
include make.cmdx.inc

SRC_DIRS := ./
PY_SRC_DIR := ./
PYTHON := $(or ${PYTHON}, python)

QEMODS := $(wildcard ${QEDIR}/*/*.a) $(wildcard ${QEDIR}/*/*/*.a)
IGNOREMODS := ${QEDIR}/LR_Modules/liblrmod.a ${QEDIR}/PW/src/libpw.a ${QEDIR}/Modules/libqemod.a
IGNOREMODS += ${QEDIR}/CPV/src/libcp.a ${QEDIR}/PHonon/Gamma/libphcg.a ${QEDIR}/GWW/bse/libbse.a ${QEDIR}/GWW/pw4gww/libpw4gww.a

QE_OBJS := $(filter-out ${QEDIR}/atomic/src/ld1.o, $(wildcard ${QEDIR}/atomic/src/*.o))
QE_OBJS += $(filter-out ${QEDIR}/PWCOND/src/condmain.o, $(wildcard ${QEDIR}/PWCOND/src/*.o))
QE_OBJS += $(filter-out $(addprefix ${QEDIR}/KCW/PP/, compute_self_hartree.o interp_r_to_k.o ups.o), $(wildcard ${QEDIR}/KCW/PP/*.o))

QEMODS := ${QEDIR}/PW/src/libpw.a ${QEDIR}/Modules/libqemod.a $(filter-out ${IGNOREMODS}, $(QEMODS))
QEMODS += ${QE_OBJS}
QEMODS += ${QEDIR}/LR_Modules/liblrmod.a ${QEDIR}/PW/src/libpw.a ${QEDIR}/Modules/libqemod.a ${QEDIR}/upflib/libupf.a
$(info 'QEMODS', ${QEMODS})

QEINC := $(wildcard ${QEDIR}/*/Makefile) $(wildcard ${QEDIR}/*/*/Makefile)
QEINC := $(addprefix -I, ${QEINC})
QEINC := $(QEINC:%/Makefile=%/) -I.

TDDFT_SOURCES = tddft_module.f90
TDDFT_FILES = ${TDDFT_SOURCES}

QEPY_ADD := qepy_sys.f90 qepy_common.f90 qepy_mod.f90 \
			scatter_mod.f90 funct.f90 potinit.f90 wfcinit.f90\
			pw_restart_new.f90 read_file_new.f90 set_para_diag.f90 setup.f90


QEPY_SOURCES := qepy_setlocal.f90 qepy_v_of_rho.f90 qepy_pw2casino_write.f90 \
              qepy_pwscf.f90 qepy_run_pwscf.f90 qepy_electrons.f90 \
              qepy_hinit1.f90 qepy_forces.f90 qepy_stop_run.f90 qepy_stress.f90 \
			  qepy_v_of_rho_all.f90 qepy_delta_e.f90

QEPY_TDDFT_SOURCES := qepy_tddft_common.f90 \
                    qepy_molecule_optical_absorption.f90 qepy_tddft_main.f90 \
                    qepy_tddft_routines.f90 qepy_tddft_setup.f90 qepy_update_ham.f90 \
                    qepy_tddft_mod.f90

QE_FILES := $(filter-out ${QEPY_ADD} ${QEPY_QEX}, ${QE_FILES})

EXT_FILES=

ifeq ($(original), yes)
   QEPY_SOURCES := ${QEPY_ADD}
   tddft = no
else
   QEPY_SOURCES += ${QEPY_ADD}
endif

QEPY_SOURCES += ${QEPY_QEX}

ifeq ($(tddft), yes)
   QEPY_SOURCES := ${QEPY_TDDFT_SOURCES} ${QEPY_SOURCES}
   QEMODS := libtddft.a ${QEMODS}
   EXT_FILES += ${TDDFT_FILES}
endif

QEPY_FILES = ${QEPY_SOURCES}

#$(info 'files', ${QEPY_FILES})
#$(info 'files', ${QE_FILES})

QEPY_OBJS = $(QEPY_SOURCES:%.f90=%.o)

WRAP_FILES = ${QE_FILES} ${QEPY_FILES} ${EXT_FILES}

F90WRAP_FILES = f90wrap_*.f90

WRAP_FPP_FILES = $(notdir $(WRAP_FILES:%.f90=%.fpp))

$(info 'F90FLAGS', ${F90FLAGS})
F2FLAGS = $(F90FLAGS) $(QEINC)
#$(info 'F2FLAGS', ${F2FLAGS})
export LDFLAGS := ${LDFLAGS}
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
   CFLAGS := $(CFLAGS) -std=c99
endif

FPP = ${F90} -E -cpp $(F2FLAGS)

PY3_DIR := $(shell python3 -m site --user-site)

default: ${GOAL}

vpath %.f90 $(SRC_DIRS)

$(filter %.o,${QEPY_OBJS}):%.o : %.f90
	$(LD) -c $(F2FLAGS) $< -o $@

${WRAP_FPP_FILES}: ${WRAP_FILES}
	@for f in ${WRAP_FILES}; do $(FPP) $$f > $$(basename $${f%.*}).fpp; done

${F90WRAP_FILES}: ${QEPY_OBJS} ${WRAP_FPP_FILES}
	@${PYTHON} -m f90wrap --f90wrap -m qepy ${WRAP_FPP_FILES} -k $(PY_SRC_DIR)/kind_map.json \
		--init-file $(PY_SRC_DIR)/__init__.py -P

.PHONY: serial mpi clean install help

serial: ${F90WRAP_FILES}
	CFLAGS='$(CFLAGS)' \
	${PYTHON} -m f90wrap --f2py-f90wrap \
		--fcompiler=$(FCOMPILER) --build-dir . \
		-c --f90exec=$(F90) --f77exec=$(F90) \
		--f90flags='${F2FLAGS}' \
		-m _qepy ${F90WRAP_FILES} $(QEPY_OBJS) \
		${QEMODS} ${QELIBS}

mpi: ${F90WRAP_FILES}
	CFLAGS='$(CFLAGS)' \
	${PYTHON} -m f90wrap --f2py-f90wrap \
		--fcompiler=$(FCOMPILER) --build-dir . \
		-c --f90exec=$(MPIF90) --f77exec=$(MPIF90) \
		--f90flags='${F2FLAGS}' \
		-m _qepy ${F90WRAP_FILES} $(QEPY_OBJS) \
		${QEMODS} ${QELIBS}

install:
	cp -r qepy _qepy*.so ${PY3_DIR}

uninstall:
	-rm -rf ${PY3_DIR}/qepy
	-rm -rf ${PY3_DIR}/_qepy*.so

clean:
	-rm -f _qepy*.so ${F90WRAP_FILES} ${WRAP_FPP_FILES}
	-rm -rf qepy
	-rm -rf src.* .libs .f2py_f2cmap
	-rm -rf *.o *.mod

rmso:
	-rm -f _qepy*.so

test: rmso mpi uninstall
	cp -r qepy _qepy*.so ${PY3_DIR}

include make.depend
