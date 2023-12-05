include ../make.qepy.inc
include ../make.depend
QESRC_DIR = upflib
PNAME = qepy_upflib
QEPY_LIBS = -L../ -lqepy_utilxlib
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I../utilxlib -I$(QEDIR)/UtilXlib

QEPY_OBJS = 

QEPY_PROG_OBJS = 

QE_OBJS = uspp.o \
uspp_param.o

LIB_OBJS_TEMP = atom.o \
atomic_number.o \
casino_pp.o \
dom.o \
dqvan2.o \
dylmr2.o \
dylmr2_gpu.o \
gen_us_dj.o \
gen_us_dy.o \
gth.o \
init_tab_atwfc.o \
init_tab_beta.o \
init_tab_qrad.o \
init_us_0.o \
init_us_1.o \
init_us_2_base.o \
init_us_2_base_gpu.o \
init_us_b0.o \
interp_atwfc.o \
interp_atwfc_gpu.o \
paw_variables.o \
pseudo_types.o \
qvan2.o \
qvan2_gpu.o \
radial_grids.o \
read_cpmd.o \
read_fhi.o \
read_ncpp.o \
read_ps.o \
read_upf_new.o \
read_upf_v1.o \
read_uspp.o \
simpsn.o \
simpsn_gpu.o \
sph_bes.o \
sph_bes_gpu.o \
sph_ind.o \
spinor.o \
splinelib.o \
upf_acc_interfaces.o \
upf_auxtools.o \
upf_const.o \
upf_error.o \
upf_invmat.o \
upf_io.o \
upf_ions.o \
upf_kinds.o \
upf_parallel_include.o \
upf_params.o \
upf_spinorb.o \
upf_to_internal.o \
upf_utils.o \
uspp.o \
uspp_data.o \
uspp_param.o \
write_upf_new.o \
wxml.o \
xmltools.o \
ylmr2.o \
ylmr2_gpu.o

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
