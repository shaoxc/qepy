include ../make.qepy.inc
include ../make.depend
QESRC_DIR = UtilXlib
PNAME = qepy_utilxlib
QEPY_LIBS = -L../ 
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) 

QEPY_OBJS = 

QEPY_PROG_OBJS = 

QEPY_FIX_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = c_mkdir.o \
clib_wrappers.o \
clocks_handler.o \
copy.o \
cptimer.o \
data_buffer.o \
device_helper.o \
divide.o \
error_handler.o \
eval_infix.o \
export_gstart_2_solvers.o \
find_free_unit.o \
fletcher32.o \
fletcher32_mod.o \
md5.o \
md5_from_file.o \
mem_counter.o \
memstat.o \
memusage.o \
mp.o \
mp_bands_util.o \
mp_base.o \
mp_base_gpu.o \
nvtx_wrapper.o \
parallel_include.o \
print_mem.o \
ptrace.o \
set_mpi_comm_4_solvers.o \
thread_util.o \
util_param.o

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
