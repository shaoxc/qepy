include ../make.qepy.inc
include ../make.depend
QESRC_DIR = KS_Solvers
PNAME = qepy_ks_solvers
QEPY_LIBS = -L../ -lqepy_utilxlib -lqepy_laxlib -lqepy_modules -lqepy_pw
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I../utilxlib -I$(QEDIR)/UtilXlib -I../laxlib -I$(QEDIR)/LAXlib -I../modules -I$(QEDIR)/Modules -I../pw -I$(QEDIR)/PW/src

QEPY_OBJS = 

QEPY_PROG_OBJS = 

QEPY_FIX_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = CG/ccgdiagg.o \
CG/ccgdiagg_gpu.o \
CG/rcgdiagg.o \
CG/rcgdiagg_gpu.o \
Davidson/cegterg.o \
Davidson/cegterg_gpu.o \
Davidson/regterg.o \
Davidson/regterg_gpu.o \
Davidson_RCI/david_rci.o \
DENSE/gram_schmidt_gamma.o \
DENSE/gram_schmidt_gamma_gpu.o \
DENSE/gram_schmidt_k.o \
DENSE/gram_schmidt_k_gpu.o \
DENSE/rotate_HSpsi_gamma.o \
DENSE/rotate_HSpsi_k.o \
DENSE/rotate_driver.o \
DENSE/rotate_driver_cuf.o \
DENSE/rotate_wfc_gamma.o \
DENSE/rotate_wfc_gamma_gpu.o \
DENSE/rotate_wfc_k.o \
DENSE/rotate_wfc_k_gpu.o \
DENSE/rotate_xpsi_gamma.o \
DENSE/rotate_xpsi_gamma_gpu.o \
DENSE/rotate_xpsi_k.o \
DENSE/rotate_xpsi_k_gpu.o \
ParO/bpcg_gamma.o \
ParO/bpcg_k.o \
ParO/paro_gamma.o \
ParO/paro_gamma_new.o \
ParO/paro_k.o \
ParO/paro_k_new.o \
ParO/pcg_gamma.o \
ParO/pcg_k.o \
PPCG/generic_cublas.o \
PPCG/ppcg_gamma.o \
PPCG/ppcg_gamma_gpu.o \
PPCG/ppcg_k.o \
PPCG/ppcg_k_gpu.o \
RMM/crmmdiagg.o \
RMM/crmmdiagg_gpu.o \
RMM/rrmmdiagg.o \
RMM/rrmmdiagg_gpu.o

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
