QESRC_DIR = LAXlib
PNAME = qepy_laxlib
QEPY_LIBS = -L$(BUILDDIR) -lqepy_utilxlib
QEPY_INCS = -I. -I$(QEDIR)/$(QESRC_DIR) -I$(BUILDDIR)/utilxlib -I$(QEDIR)/UtilXlib

QEPY_OBJS = 

QEPY_PROG_OBJS = 

QE_OBJS = mp_diag.o \
la_helper.o

LIB_OBJS_TEMP = cdiaghg.o \
distools.o \
dspev_drv.o \
la_error.o \
la_helper.o \
la_module.o \
la_param.o \
la_types.o \
mp_diag.o \
ptoolkit.o \
rdiaghg.o \
transto.o \
zhpev_drv.o

include ../make.qepy.inc
include ../make.depend
