include ../make.qepy.inc
include ../make.depend
QESRC_DIR = FFTXlib/src
PNAME = qepy_fftxlib
QEPY_LIBS = -L../ -lqepy_utilxlib
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I../utilxlib -I$(QEDIR)/UtilXlib

QEPY_OBJS = 

QEPY_PROG_OBJS = 

QE_OBJS = fft_types.o

LIB_OBJS_TEMP = fft_buffers.o \
fft_error.o \
fft_fwinv.o \
fft_ggen.o \
fft_helper_subroutines.o \
fft_interfaces.o \
fft_interpolate.o \
fft_parallel.o \
fft_parallel_2d.o \
fft_param.o \
fft_scalar.DFTI.o \
fft_scalar.ESSL.o \
fft_scalar.FFTW.o \
fft_scalar.FFTW3.o \
fft_scalar.SX6.o \
fft_scalar.cuFFT.o \
fft_scalar.o \
fft_scatter.o \
fft_scatter_2d.o \
fft_scatter_2d_gpu.o \
fft_scatter_gpu.o \
fft_smallbox.o \
fft_smallbox_type.o \
fft_stick.o \
fft_support.o \
fft_types.o \
fftw.o \
fftw_dp.o \
fftw_interfaces.o \
fftw_sp.o \
scatter_mod.o \
stick_base.o \
tg_gather.o

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
