include ../make.qepy.inc
include ../make.depend
QESRC_DIR = XClib
PNAME = qepy_xclib
QEPY_LIBS = -L../ -lqepy_upflib
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I../upflib -I$(QEDIR)/upflib

QEPY_OBJS = 

QEPY_PROG_OBJS = 

QEPY_FIX_OBJS = 

QE_OBJS = dft_setting_routines.o

LIB_OBJS_TEMP = beefun.o \
dft_setting_params.o \
dft_setting_routines.o \
pbecor.o \
qe_constants.o \
qe_dft_list.o \
qe_drivers_d_gga.o \
qe_drivers_d_lda_lsda.o \
qe_drivers_gga.o \
qe_drivers_lda_lsda.o \
qe_drivers_mgga.o \
qe_funct_corr_gga.o \
qe_funct_corr_lda_lsda.o \
qe_funct_exch_gga.o \
qe_funct_exch_lda_lsda.o \
qe_funct_mgga.o \
qe_kind.o \
xc_beef_interface.o \
xc_lib.o \
xc_wrapper_d_gga.o \
xc_wrapper_d_lda_lsda.o \
xc_wrapper_gga.o \
xc_wrapper_lda_lsda.o \
xc_wrapper_mgga.o \
xclib_error.o \
xclib_utils_and_para.o

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
