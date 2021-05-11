QEDIR := $(or ${qedir}, ../../../)

include ${QEDIR}/make.inc

SRC_DIRS := ./
PY_SRC_DIR := ./

QEMODS= ${QEDIR}/PW/src/libpw.a ${QEDIR}/Modules/libqemod.a ${QEDIR}/KS_Solvers/libks_solvers.a ${QEDIR}/FFTXlib/libqefft.a \
	   ${QEDIR}/LAXlib/libqela.a ${QEDIR}/UtilXlib/libutil.a ${QEDIR}/dft-d3/libdftd3qe.a

QEINC =-I${QEDIR}/Modules/ -I${QEDIR}/KS_Solvers/ -I${QEDIR}/FFTXlib/ \
	   -I${QEDIR}/LAXlib/ -I${QEDIR}/UtilXlib/ -I${QEDIR}/dft-d3/ -I${QEDIR}/PW/src/ -I. -I${QEDIR}/iotk/src/

MODULES_SOURCES = constants.f90 cell_base.f90 ions_base.f90
MODULES_FILES = $(addprefix ${QEDIR}/Modules/,${MODULES_SOURCES})

PW_SOURCES = pwcom.f90 scf_mod.f90 read_file_new.f90 punch.f90 atomic_wfc_mod.f90
PW_FILES = $(addprefix ${QEDIR}/PW/src/,${PW_SOURCES})

QEPY_SOURCES := qepy_scatter_mod.f90 \
              qepy_common.f90 qepy_mod.f90 \
              qepy_setlocal.f90 qepy_v_of_rho.f90 qepy_pw2casino_write.f90 \
              qepy_hinit1.f90 qepy_potinit.f90 qepy_wfcinit.f90 qepy_pw_restart_new.f90 \
              qepy_init_run.f90 qepy_pwscf.f90 qepy_run_pwscf.f90 qepy_electrons.f90 \
              qepy_forces.f90 qepy_stop_run.f90

OLDXML_SOURCES := oldxml_qexml.f90 oldxml_xml_io_base.f90 \
                oldxml_io_rho_xml.f90 oldxml_pw_restart.f90 \
                oldxml_wfcinit.f90 oldxml_potinit.f90 oldxml_read_file.f90

FOLDXML =
ifeq ($(oldxml), yes)
   QEPY_SOURCES := ${OLDXML_SOURCES} ${QEPY_SOURCES}
   FOLDXML = -D__OLDXML
endif
QEPY_FILES = ${QEPY_SOURCES}

#$(info 'files', ${QEPY_FILES})

QEPY_OBJS = $(QEPY_SOURCES:%.f90=%.o)

WRAP_FILES = ${MODULES_FILES} ${PW_FILES} ${QEPY_FILES}

F90WRAP_FILES = f90wrap_*.f90

WRAP_FPP_FILES = $(notdir $(WRAP_FILES:%.f90=%.fpp))

PWFLAGS = $(F90FLAGS) $(QEINC) $(FOLDXML)
NOTI = -fPIC -nomodule -qopenmp -fpp -mcmodel=large -g -C
F2FLAGS = $(filter-out $(NOTI),$(PWFLAGS)) $(FOLDXML)

$(info 'F2FLAGS', ${F2FLAGS})

ifeq ($(findstring -D__MPI, ${F2FLAGS}), -D__MPI)
   GOAL=mpi
else
   GOAL=serial
endif


FPP = ${F90} -E -cpp $(F2FLAGS)

ifeq ($(prefix), )
    PY2_DIR := $(shell python -m site --user-site)
    PY3_DIR := $(shell python3 -m site --user-site)
ifeq ($(PY3_DIR), )
    PY3_DIR = ${PY2_DIR}
else
    PY3_DIR = $(shell python3 -m site --user-site)
endif
else
    PY2_DIR = $(prefix)
    PY3_DIR = $(prefix)
endif

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
	f2py-f90wrap --fcompiler=intelem --build-dir . \
		--opt=-O2 \
		-c -m _qepy ${F90WRAP_FILES} $(QEPY_OBJS) \
		${F2FLAGS} ${QEMODS} $(LIBOBJS) $(QELIBS) -liomp5

mpi: ${F90WRAP_FILES}
	f2py-f90wrap \
		--fcompiler=intelem --build-dir . \
		--compiler=intelem \
		-c --f90exec=mpiifort --f77exec=mpiifort \
		--opt=-O2 \
		-m _qepy ${F90WRAP_FILES} $(QEPY_OBJS) \
		${F2FLAGS} ${QEMODS} $(LIBOBJS) $(QELIBS)


install:
	cp -r qepy _qepy*.so ${PY2_DIR}
	cp -r qepy _qepy*.so ${PY3_DIR}

uninstall:
	-rm -rf ${PY2_DIR}/qepy
	-rm -rf ${PY3_DIR}/qepy
	-rm -rf ${PY2_DIR}/_qepy*.so
	-rm -rf ${PY3_DIR}/_qepy*.so

clean:
	-rm -f _qepy*.so ${F90WRAP_FILES} ${WRAP_FPP_FILES}
	-rm -rf qepy
	-rm -rf src.* .libs .f2py_f2cmap
	-rm -rf *.o *.mod

help:
	@echo ""
	@echo "make help    : show this help information"
	@echo "make serial  : build serial version (default)"
	@echo "make mpi     : build parallel [MPI] version"
	@echo "make install : install to given directory (default is user site-packages)"
	@echo "make clean   : remove the package form given directory (default is user site-packages)"
	@echo ""
	@echo "Variables:"
	@echo "prefix       : given the folder for installation"
	@echo "oldxml       : support read old version xml with iotk"

include make.depend
