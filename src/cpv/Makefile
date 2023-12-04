include ../make.qepy.inc
include ../make.depend
QESRC_DIR = CPV/src
PNAME = qepy_cpv
QEPY_LIBS = -L../ $(MODS_BASE)
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) $(INCS_MODS_BASE) $(INCS_MODS_BASE_QE)

QEPY_OBJS = 

QEPY_PROG_OBJS = cppp.o \
cprstart.o \
manycp.o \
wfdd.o

QEPY_FIX_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = berryion.o \
bforceion.o \
cell_nose.o \
cg.o \
cg_sub.o \
chargedensity.o \
cp_autopilot.o \
cp_emass.o \
cp_interfaces.o \
cp_restart_new.o \
cplib.o \
cplib_meta.o \
cpr.o \
cpr_loop.o \
cpr_mod.o \
dealloc.o \
dforceb.o \
efermi.o \
efield.o \
eigs0.o \
electrons.o \
electrons_nose.o \
energies.o \
ensemble_dft.o \
entropy.o \
environ_cp_module.o \
exch_corr.o \
exx_cg.o \
exx_es.o \
exx_gs.o \
exx_module.o \
exx_pair.o \
exx_psi.o \
exx_vofr.o \
forces.o \
fromscra.o \
gram.o \
gtable.o \
init.o \
init_run.o \
inner_loop_cold.o \
input.o \
ions_nose.o \
ions_positions.o \
ksstates.o \
ldaU.o \
ldaU_cp.o \
ldaUpen.o \
mainvar.o \
makov_payne.o \
metaxc.o \
modules.o \
move_electrons.o \
newd.o \
nl_base.o \
nlcc.o \
ortho.o \
ortho_base.o \
phasefactor.o \
plugin_add_potential.o \
plugin_clean.o \
plugin_clock.o \
plugin_energy.o \
plugin_ext_forces.o \
plugin_get_potential.o \
plugin_init_base.o \
plugin_init_cell.o \
plugin_init_ions.o \
plugin_initialization.o \
plugin_int_forces.o \
plugin_print_energies.o \
plugin_print_info.o \
plugin_read_input.o \
plugin_utilities.o \
potentials.o \
pres_ai_mod.o \
print_out.o \
printout_base.o \
problem_size.o \
pseudo_base.o \
pseudopot.o \
pseudopot_sub.o \
qmatrixd.o \
qqberry.o \
restart.o \
restart_sub.o \
runcp.o \
sic.o \
smallbox.o \
smallbox_grid.o \
smallbox_gvec.o \
smallbox_lib.o \
smallbox_subs.o \
spinsq.o \
spline.o \
stop_run.o \
stress.o \
vofrho.o \
vol_clu.o \
wannier.o \
wannier_base.o \
wave.o \
wave_base.o \
wave_types.o \
wf.o

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
