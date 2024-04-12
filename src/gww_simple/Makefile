QESRC_DIR = GWW/simple
PNAME = qepy_gww_simple
QEPY_LIBS = -L$(BUILDDIR) -lqepy_gww_gww -lqepy_gww_pw4gww $(MODS_PW)
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I$(BUILDDIR)/gww_gww -I$(QEDIR)/GWW/gww -I$(BUILDDIR)/gww_pw4gww -I$(QEDIR)/GWW/pw4gww $(INCS_MODS_PW) $(INCS_MODS_PW_QE)

QEPY_OBJS = 

QEPY_PROG_OBJS = simple.o

QE_OBJS = 

LIB_OBJS_TEMP = commutator.o \
epe.o \
gk_sort_limit.o \
init_us_2_max.o \
input_simple.o \
khamiltonian.o \
openfile_simple.o \
product_basis.o \
read_export.o \
stop_pp.o \
v_product.o \
wfc_basis.o

include ../make.qepy.inc
include ../make.depend
