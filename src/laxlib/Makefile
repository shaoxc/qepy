include ../make.qepy.inc
include ../make.depend
QESRC_DIR = LAXlib
PNAME = qepy_laxlib
QEPY_LIBS = -L../ -lqepy_utilxlib
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I../utilxlib -I$(QEDIR)/UtilXlib

QEPY_OBJS = 

QEPY_PROG_OBJS = 

QE_OBJS = mp_diag.o \
la_helper.o

LIB_OBJS_TEMP = cdiaghg.o \
distools.o \
dspev_drv.o \
la_error.o \
la_helper.o \
la_module.o \
la_param.o \
la_types.o \
mp_diag.o \
ptoolkit.o \
rdiaghg.o \
transto.o \
zhpev_drv.o

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
