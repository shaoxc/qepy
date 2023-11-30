include ../make.qepy.inc
include ../make.depend
QESRC_DIR = GWW/simple_bse
PNAME = qepy_gww_simple_bse
QEPY_LIBS = -L../ -lqepy_gww_gww $(MODS_BASE)
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I../gww_gww -I$(QEDIR)/GWW/gww $(INCS_MODS_BASE) $(INCS_MODS_BASE_QE)

QEPY_OBJS = 

QEPY_PROG_OBJS = simple_bse.o

QEPY_FIX_OBJS = 

QE_OBJS = 

LIB_OBJS_TEMP = build_eemat.o \
derived_object.o \
diago_exc_cg.o \
diago_exc_sd.o \
hamiltonian.o \
input_simple_exc.o \
lanczos.o \
simple_eigen.o \
simple_object.o \
spectrum.o \
start_end.o

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
