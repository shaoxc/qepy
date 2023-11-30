include ../make.qepy.inc
include ../make.depend
QESRC_DIR = HP/src
PNAME = qepy_hp
QEPY_LIBS = -L../ $(MODS_LR)
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) $(INCS_MODS_LR) $(INCS_MODS_LR_QE)

QEPY_OBJS = 

QEPY_PROG_OBJS = hp_main.o

QEPY_FIX_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = hp_R_points.o \
hp_allocate_q.o \
hp_bcast_input.o \
hp_calc_chi.o \
hp_check_pert.o \
hp_check_type.o \
hp_clean_q.o \
hp_close_q.o \
hp_dealloc_1.o \
hp_dealloc_2.o \
hp_dealloc_q.o \
hp_dnsq.o \
hp_dnstot_sum_q.o \
hp_dvpsi_pert.o \
hp_find_inequiv_sites.o \
hp_generate_grids.o \
hp_init.o \
hp_init_q.o \
hp_load_q.o \
hp_ns_trace.o \
hp_openfil_q.o \
hp_postproc.o \
hp_prepare_q.o \
hp_print_clock.o \
hp_psymdvscf.o \
hp_q_points.o \
hp_read_chi.o \
hp_read_dnsq.o \
hp_readin.o \
hp_rotate_dnsq.o \
hp_run_nscf.o \
hp_setup_q.o \
hp_solve_linear_system.o \
hp_stop_smoothly.o \
hp_summary.o \
hp_summary_q.o \
hp_symdnsq.o \
hp_symdvscf.o \
hp_write_chi.o \
hp_write_chi_full.o \
hp_write_dnsq.o \
hpcom.o

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
