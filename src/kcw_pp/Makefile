include ../make.qepy.inc
include ../make.depend
QESRC_DIR = KCW/PP
PNAME = qepy_kcw_pp
QEPY_LIBS = -L../ -lqepy_kcw $(MODS_LR)
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I../kcw -I$(QEDIR)/KCW/src $(INCS_MODS_LR) $(INCS_MODS_LR_QE)

QEPY_OBJS = 

QEPY_PROG_OBJS = compute_self_hartree.o \
interp_r_to_k.o

QE_OBJS = 

LIB_OBJS_TEMP = bcast_kcw_pp_input.o \
input_pp_summary.o \
kcw_bands.o \
kcw_pp_readin.o \
read_hr.o \
sh_setup.o

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
