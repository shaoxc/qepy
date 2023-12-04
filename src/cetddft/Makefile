include ../make.qepy.inc
include ../make.depend
QESRC_DIR = ../ce-tddft/src
PNAME = qepy_cetddft
QEPY_LIBS = -L../ $(MODS_LR)
QEPY_INCS = -I. $(INCS_MODS_LR) $(INCS_MODS_LR_QE)

QEPY_OBJS = qepy_molecule_optical_absorption.o \
qepy_tddft_common.o \
qepy_tddft_main.o \
qepy_tddft_mod.o \
qepy_tddft_routines.o \
qepy_tddft_setup.o \
qepy_update_ham.o

QEPY_PROG_OBJS = tddft_main.o

QEPY_FIX_OBJS = tddft_version.o tddft_module.o tddft_routines.o \
	     tddft_setup.o tddft_ch_psi_all.o tddft_cgsolver.o \
	     update_ham.o apply_efield.o molecule_optical_absorption.o \
	     molecule_operators.o stop_code.o trajectory.o wavepacket.o

QE_OBJS = 

LIB_OBJS_TEMP =

vpath %.f90 ../fix/ ../cmdx/ $(SRC_DIR) $(QESRC_DIR)

QEPY_OBJS := ${QEPY_OBJS} ${QEPY_FIX_OBJS} ${QEPY_PROG_OBJS}
$(info $(QEPY_OBJS))
LIB_OBJS_TEMP := $(filter-out ${QEPY_OBJS}, ${LIB_OBJS_TEMP})
LIB_OBJS := $(addprefix $(QESRC_DIR)/,$(LIB_OBJS_TEMP))

ifeq ($(strip $(WRAP_OBJS)),)
QEPY_OBJS = qepy_fake.o
endif

.PHONY: all

${F90WRAP_FILES}: ${QEPY_OBJS} ${WRAP_FPP_FILES}
	@${PYTHON} -m f90wrap --f90wrap -m ${PNAME} --f90-mod-name lib${PNAME} ${WRAP_FPP_FILES} -k $(PY_SRC_DIR)/kind_map.json -P 
