include ../make.qepy.inc
include ../make.depend
QESRC_DIR = LR_Modules
PNAME = qepy_lr_modules
QEPY_LIBS = -L../ $(MODS_PW)
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) $(INCS_MODS_PW) $(INCS_MODS_PW_QE)

QEPY_OBJS = 

QEPY_PROG_OBJS = 

QEPY_FIX_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = Coul_cut_2D_ph.o \
adddvepsi_us.o \
adddvhubscf.o \
adddvscf.o \
addusdbec.o \
addusdbec_nc.o \
addusddenseq.o \
addusldos.o \
apply_dpot.o \
apply_dpot_mod.o \
ccg_psi.o \
ccgsolve_all.o \
cft_wave.o \
cg_psi.o \
cgsolve_all.o \
ch_psi_all.o \
ch_psi_all_complex.o \
check_q_points_sym.o \
compute_vsgga.o \
dfpt_tetra_mod.o \
dgradcorr.o \
dnonloccorr.o \
dv_of_drho.o \
dv_rVV10.o \
dv_vdW_DF.o \
dynmat_sub.o \
efermi_shift.o \
h_prec.o \
incdrhoscf.o \
incdrhoscf_nc.o \
lanczos_nonhermitian.o \
lanczos_nonhermitian_c.o \
lanczos_pseudohermitian.o \
lanczos_pseudohermitian_c.o \
localdos.o \
lr_addusddens.o \
lr_dot.o \
lr_dot_magnons.o \
lr_orthoUwfc.o \
lr_sm1_psi.o \
lr_sym_mod.o \
lrcom.o \
mix_pot.o \
mode_group.o \
newdq.o \
orthogonalize.o \
qdipol_cryst.o \
response_kernels.o \
set_dbecsum_nc.o \
set_int3_nc.o \
set_intq_nc.o \
set_kplusq_nc.o \
set_small_group_of_q.o \
setqmod.o \
setup_alpha_pv.o \
setup_dgc.o \
setup_dmuxc.o \
setup_nbnd_occ.o \
setup_nscf.o \
sgam_lr.o \
smallgq.o \
star_q.o

vpath %.f90 ../fix/ ../cmdx/ $(SRC_DIR) $(QEDIR)/$(QESRC_DIR)

QEPY_OBJS := ${QEPY_OBJS} ${QEPY_FIX_OBJS} ${QEPY_PROG_OBJS}
LIB_OBJS_TEMP := $(filter-out ${QEPY_OBJS}, ${LIB_OBJS_TEMP})
LIB_OBJS := $(addprefix $(QEDIR)/$(QESRC_DIR)/,$(LIB_OBJS_TEMP))

ifeq ($(strip $(WRAP_OBJS)),)
QEPY_OBJS = qepy_fake.o
endif

.PHONY: all

all : ${F90WRAP_FILES}
	${F2PY} -m lib${PNAME} ${F90WRAP_FILES} $(QEPY_OBJS) $(LIB_OBJS) $(QEPY_LIBS) ${QELIBS}
	cp lib$(PNAME).*.so ../lib$(PNAME).so
	if test -d $(PNAME); then mv $(PNAME) ../; fi

${F90WRAP_FILES}: ${QEPY_OBJS} ${WRAP_FPP_FILES}
	@${PYTHON} -m f90wrap --f90wrap -m ${PNAME} --f90-mod-name lib${PNAME} ${WRAP_FPP_FILES} -k $(PY_SRC_DIR)/kind_map.json -P 
