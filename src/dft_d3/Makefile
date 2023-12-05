include ../make.qepy.inc
include ../make.depend
QESRC_DIR = dft-d3
PNAME = qepy_dft_d3
QEPY_LIBS = -L../ -lqepy_modules -lqepy_utilxlib -lqepy_laxlib
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I../modules -I$(QEDIR)/Modules -I../utilxlib -I$(QEDIR)/UtilXlib -I../laxlib -I$(QEDIR)/LAXlib

QEPY_OBJS = 

QEPY_PROG_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = api.o \
common.o \
core.o \
dftd3_qe.o \
pars.o \
sizes.o

vpath %.f90 ../fix/ ../cmdx/ $(SRC_DIR) $(QEDIR)/$(QESRC_DIR)

QEPY_OBJS := ${QEPY_OBJS} ${QEPY_FIX_OBJS} ${QEPY_PROG_OBJS}
LIB_OBJS_TEMP := $(filter-out ${QEPY_OBJS}, ${LIB_OBJS_TEMP})
LIB_OBJS := $(addprefix $(QEDIR)/$(QESRC_DIR)/,$(LIB_OBJS_TEMP))

ifeq ($(strip $(WRAP_OBJS)),)
QEPY_OBJS = qepy_fake.o
endif

.PHONY: all

${F90WRAP_FILES}: ${QEPY_OBJS} ${WRAP_FPP_FILES}
	@${PYTHON} -m f90wrap --f90wrap -m ${PNAME} --f90-mod-name lib${PNAME} ${WRAP_FPP_FILES} ${CROSS_FPP_FILES} -k $(PY_SRC_DIR)/kind_map.json -P 
