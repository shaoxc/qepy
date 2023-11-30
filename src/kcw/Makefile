include ../make.qepy.inc
include ../make.depend
QESRC_DIR = KCW/src
PNAME = qepy_kcw
QEPY_LIBS = -L../ $(MODS_LR)
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) $(INCS_MODS_LR) $(INCS_MODS_LR_QE)

QEPY_OBJS = 

QEPY_PROG_OBJS = kcw.o

QEPY_FIX_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = alpha_corr.o \
apply_u_matrix.o \
bare_pot.o \
bcast_kcw_input.o \
bcast_wfc.o \
calculate_phase.o \
clean_pw_kcw.o \
close_kcw.o \
compute_map_ikq_single.o \
convert_kpts_names.o \
coulomb.o \
full_ham.o \
group_orbitals.o \
ham_R0_2nd.o \
input_summary.o \
interpolation.o \
kcw_R_points.o \
kcw_allocate_q.o \
kcw_comm.o \
kcw_deallocate_q.o \
kcw_ham.o \
kcw_init_q.o \
kcw_initialize_ph.o \
kcw_io_new.o \
kcw_openfilq.o \
kcw_prepare_q.o \
kcw_q_setup.o \
kcw_readin.o \
kcw_run_nscf.o \
kcw_screen.o \
kcw_setup.o \
kcw_setup_ham.o \
kcw_setup_screen.o \
koopmans_ham.o \
ks_hamiltonian.o \
print_clock_kcw.o \
read_alpha.o \
read_wannier.o \
rho_of_q.o \
rotate_ks.o \
screen_coeff.o \
self_hartree.o \
setup_coulomb_exx.o \
solve_linter_koop_mod.o \
structure_factor.o \
wann2kcw.o \
write_hr_to_file.o

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
