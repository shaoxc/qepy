include ../make.qepy.inc
include ../make.depend
QESRC_DIR = GWW/head
PNAME = qepy_gww_head
QEPY_LIBS = -L../ -lqepy_gww_minpack -lqepy_phonon_ph $(MODS_LR)
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I../gww_minpack -I$(QEDIR)/GWW/minpack -I../phonon_ph -I$(QEDIR)/PHonon/PH $(INCS_MODS_LR) $(INCS_MODS_LR_QE)

QEPY_OBJS = 

QEPY_PROG_OBJS = head.o

QEPY_FIX_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = bcast_gwl_input.o \
close_gwl.o \
gwl_readin.o \
lanczos_k.o \
openfilq_gwl.o \
solve_head.o

vpath %.f90 ../fix/ ../cmdx/ $(SRC_DIR) $(QEDIR)/$(QESRC_DIR)

QEPY_OBJS := ${QEPY_OBJS} ${QEPY_FIX_OBJS} ${QEPY_PROG_OBJS}
LIB_OBJS_TEMP := $(filter-out ${QEPY_OBJS}, ${LIB_OBJS_TEMP})
LIB_OBJS := $(addprefix $(QEDIR)/$(QESRC_DIR)/,$(LIB_OBJS_TEMP))

ifeq ($(strip $(WRAP_OBJS)),)
QEPY_OBJS = qepy_fake.o
endif

.PHONY: all

${F90WRAP_FILES}: ${QEPY_OBJS} ${WRAP_FPP_FILES}
	@${PYTHON} -m f90wrap --f90wrap -m ${PNAME} --f90-mod-name lib${PNAME} ${WRAP_FPP_FILES} -k $(PY_SRC_DIR)/kind_map.json -P 
