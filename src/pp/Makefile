include ../make.qepy.inc
include ../make.depend
QESRC_DIR = PP/src
PNAME = qepy_pp
QEPY_LIBS = -L../ $(MODS_PW)
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) $(INCS_MODS_PW) $(INCS_MODS_PW_QE)

QEPY_OBJS = 

QEPY_PROG_OBJS = average.o \
band_interpolation.o \
bands.o \
d3hess.o \
dos.o \
epsilon.o \
fermi_proj.o \
fermi_velocity.o \
fermisurface.o \
initial_state.o \
molecularpdos.o \
open_grid.o \
oscdft_et.o \
oscdft_pp.o \
plan_avg.o \
plotband.o \
plotproj.o \
plotrho.o \
poormanwannier.o \
postproc.o \
postrism.o \
ppacf.o \
projwfc.o \
pw2bgw.o \
pw2critic.o \
pw2gw.o \
pw2wannier90.o \
sumpdos.o \
wannier_ham.o \
wannier_plot.o \
wfck2r.o

QEPY_FIX_OBJS = 

QE_OBJS = elf.o \
local_dos.o \
local_dos_mag.o

LIB_OBJS_TEMP = add_shift_cc.o \
add_shift_lc.o \
add_shift_us.o \
adduscore.o \
addusdens1d.o \
atomic_wfc_nc_proj.o \
chdens_bspline.o \
chdens_bspline_sphere.o \
chdens_module.o \
compute_ppsi.o \
compute_sigma_avg.o \
cube.o \
d_matrix_nc.o \
d_matrix_so.o \
do_initial_state.o \
do_shift_ew.o \
dosg.o \
elf.o \
fermisurfer_common.o \
fouriermod.o \
ggen1d.o \
globalmod.o \
idwmod.o \
local_dos.o \
local_dos1d.o \
local_dos_mag.o \
openfil_pp.o \
oscdft_et_mod.o \
oscdft_pp_mod.o \
partialdos.o \
paw_postproc.o \
projections_mod.o \
projwfc_box.o \
punch_plot.o \
punch_rism.o \
read_proj.o \
smallgk.o \
solvdens.o \
stm.o \
stop_pp.o \
sum_band_kin.o \
sym_band.o \
vasp_read_chgcar_mod.o \
vasp_xml_module.o \
wannier_enrg.o \
wannier_proj.o \
work_function.o \
write_hamiltonians.o \
write_io_header.o \
write_p_avg.o \
write_proj.o \
xc_vdW_scale_mod.o

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
