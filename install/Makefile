QEDIR := $(or ${qedir}, ../../../)

include ${QEDIR}/make.inc

SRC_DIRS := ./
PY_SRC_DIR := ./

QEMODS := $(wildcard ${QEDIR}/*/*.a)
QEMODS := ${QEDIR}/PW/src/libpw.a ${QEDIR}/Modules/libqemod.a $(filter-out ${QEDIR}/Modules/libqemod.a, $(QEMODS))
#QEMODS := ${QEDIR}/PW/src/libpw.a ${QEDIR}/Modules/libqemod.a ${QEDIR}/KS_Solvers/libks_solvers.a ${QEDIR}/FFTXlib/libqefft.a \
#       ${QEDIR}/LAXlib/libqela.a ${QEDIR}/UtilXlib/libutil.a ${QEDIR}/dft-d3/libdftd3qe.a
#
#$(info 'qemods', ${QEMODS})

QEINC := $(wildcard ${QEDIR}/*/Makefile) $(wildcard ${QEDIR}/*/src/Makefile)
QEINC := $(addprefix -I, ${QEINC})
QEINC := $(QEINC:%/Makefile=%/) -I.

MODULES_SOURCES = constants.f90 cell_base.f90 ions_base.f90 wavefunctions.f90 funct.f90 \
				  recvec.f90 control_flags.f90
MODULES_FILES = $(addprefix ${QEDIR}/Modules/,${MODULES_SOURCES})

PW_SOURCES = pwcom.f90 scf_mod.f90 read_file_new.f90 punch.f90 \
			 atomic_wfc_mod.f90 close_files.f90 stress.f90
PW_FILES = $(addprefix ${QEDIR}/PW/src/,${PW_SOURCES})

TDDFT_SOURCES = tddft_module.f90
TDDFT_FILES = ${TDDFT_SOURCES}

QEPY_ADD := qepy_scatter_mod.f90 qepy_common.f90 qepy_mod.f90

QEPY_SOURCES := qepy_setlocal.f90 qepy_v_of_rho.f90 qepy_pw2casino_write.f90 \
              qepy_hinit1.f90 qepy_potinit.f90 qepy_wfcinit.f90 qepy_pw_restart_new.f90 \
              qepy_init_run.f90 qepy_pwscf.f90 qepy_run_pwscf.f90 qepy_electrons.f90 \
              qepy_forces.f90 qepy_stop_run.f90 qepy_api.f90

OLDXML_SOURCES := oldxml_qexml.f90 oldxml_xml_io_base.f90 \
                oldxml_io_rho_xml.f90 oldxml_pw_restart.f90 \
                oldxml_wfcinit.f90 oldxml_potinit.f90 oldxml_read_file.f90

QEPY_TDDFT_SOURCES := qepy_tddft_common.f90 \
                    qepy_molecule_optical_absorption.f90 qepy_tddft_main.f90 \
                    qepy_tddft_routines.f90 qepy_tddft_setup.f90 qepy_update_ham.f90 \
                    qepy_tddft_mod.f90

FOLDXML =
EXT_FILES=

ifeq ($(original), yes)
   QEPY_SOURCES := ${QEPY_ADD}
   oldxml = no
   tddft = no
else
   QEPY_SOURCES += ${QEPY_ADD}
endif

ifeq ($(oldxml), yes)
   QEPY_SOURCES := ${OLDXML_SOURCES} ${QEPY_SOURCES}
   FOLDXML = -D__OLDXML
endif

ifeq ($(tddft), yes)
   QEPY_SOURCES := ${QEPY_TDDFT_SOURCES} ${QEPY_SOURCES}
   QEMODS := libtddft.a ${QEMODS}
   EXT_FILES += ${TDDFT_FILES}
endif

QEPY_FILES = ${QEPY_SOURCES}

#$(info 'files', ${QEPY_FILES})

QEPY_OBJS = $(QEPY_SOURCES:%.f90=%.o)

WRAP_FILES = ${MODULES_FILES} ${PW_FILES} ${QEPY_FILES} ${EXT_FILES}

F90WRAP_FILES = f90wrap_*.f90

WRAP_FPP_FILES = $(notdir $(WRAP_FILES:%.f90=%.fpp))

PWFLAGS = $(F90FLAGS) $(QEINC) $(FOLDXML)
NOTI = 
F2FLAGS = $(filter-out $(NOTI),$(PWFLAGS))

$(info 'F2FLAGS', ${F2FLAGS})

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


FPP = ${F90} -E -cpp $(F2FLAGS)

PY3_DIR := $(shell python3 -m site --user-site)

default: ${GOAL}

vpath %.f90 $(SRC_DIRS)

$(filter %.o,${QEPY_OBJS}):%.o : %.f90
	$(LD) -c $(PWFLAGS) $< -o $@

${WRAP_FPP_FILES}: ${WRAP_FILES}
	for f in ${WRAP_FILES}; do $(FPP) $$f > $$(basename $${f%.*}).fpp; done

${F90WRAP_FILES}: ${QEPY_OBJS} ${WRAP_FPP_FILES}
	f90wrap -v -m qepy ${WRAP_FPP_FILES} -k $(PY_SRC_DIR)/kind_map \
		--init-file $(PY_SRC_DIR)/init.py -P

.PHONY: serial mpi clean install help

serial: ${F90WRAP_FILES}
	f2py-f90wrap \
		--fcompiler=$(FCOMPILER) --build-dir . \
		-c --f90exec=$(F90) --f77exec=$(F90) \
		--opt=-O2 \
		--f90flags='${F2FLAGS}' \
		-m _qepy ${F90WRAP_FILES} $(QEPY_OBJS) \
		${QEMODS} $(LIBOBJS) $(QELIBS) -liomp5

mpi: ${F90WRAP_FILES}
	f2py-f90wrap \
		--fcompiler=$(FCOMPILER) --build-dir . \
		-c --f90exec=$(MPIF90) --f77exec=$(MPIF90) \
		--f90flags='${F2FLAGS}' \
		-m _qepy ${F90WRAP_FILES} $(QEPY_OBJS) \
		${QEMODS} $(LIBOBJS) $(QELIBS)

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
