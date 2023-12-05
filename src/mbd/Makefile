include ../make.qepy.inc
include ../make.depend
QESRC_DIR = external/mbd/src
PNAME = qepy_mbd
QEPY_LIBS = -L../ 
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) 

QEPY_OBJS = 

QEPY_PROG_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = mbd.o \
mbd_c_api.o \
mbd_constants.o \
mbd_coulomb.o \
mbd_damping.o \
mbd_dipole.o \
mbd_formulas.o \
mbd_geom.o \
mbd_gradients.o \
mbd_hamiltonian.o \
mbd_lapack.o \
mbd_linalg.o \
mbd_matrix.o \
mbd_methods.o \
mbd_rpa.o \
mbd_scs.o \
mbd_ts.o \
mbd_utils.o \
mbd_vdw_param.o

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
