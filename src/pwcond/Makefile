include ../make.qepy.inc
include ../make.depend
QESRC_DIR = PWCOND/src
PNAME = qepy_pwcond
QEPY_LIBS = -L../ $(MODS_PW)
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) $(INCS_MODS_PW) $(INCS_MODS_PW_QE)

QEPY_OBJS = 

QEPY_PROG_OBJS = condmain.o

QEPY_FIX_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = allocate_cond.o \
bessj.o \
compbs.o \
compbs_2.o \
cond_out.o \
cond_restart.o \
condcom.o \
do_cond.o \
eigenchnl.o \
form_zk.o \
four.o \
free_mem.o \
gep_x.o \
gramsh.o \
hev_ab.o \
init_cond.o \
init_gper.o \
init_orbitals.o \
integrals.o \
jbloch.o \
kbloch.o \
local.o \
local_set.o \
openfil_cond.o \
plus_u_setup.o \
poten.o \
print_clock_pwcond.o \
realus_scatt.o \
rotproc.o \
save_cond.o \
scat_states_plot.o \
scatter_forw.o \
summary_band.o \
summary_tran.o \
sunitary.o \
transmit.o

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
