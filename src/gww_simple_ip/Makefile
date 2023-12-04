include ../make.qepy.inc
include ../make.depend
QESRC_DIR = GWW/simple_ip
PNAME = qepy_gww_simple_ip
QEPY_LIBS = -L../ -lqepy_gww_pw4gww $(MODS_PW)
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I../gww_pw4gww -I$(QEDIR)/GWW/pw4gww $(INCS_MODS_PW) $(INCS_MODS_PW_QE)

QEPY_OBJS = 

QEPY_PROG_OBJS = simple_ip.o

QEPY_FIX_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = diagonalization.o \
dielectric.o \
input_simple_ip.o \
interpolation.o \
simple_ip_objects.o \
start_end.o \
tetra_mod1.o

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
