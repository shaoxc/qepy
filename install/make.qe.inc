QEDIR = $(or ${qedir}, ../../../)

include ${QEDIR}/make.inc

MODULES_SOURCES = constants.f90 cell_base.f90 ions_base.f90 wavefunctions.f90 funct.f90 \
				  recvec.f90 control_flags.f90
MODULES_FILES = $(addprefix ${QEDIR}/Modules/,${MODULES_SOURCES})

PW_SOURCES = pwcom.f90 scf_mod.f90 read_file_new.f90 punch.f90 \
			 atomic_wfc_mod.f90 close_files.f90 stress.f90 electrons.f90
PW_FILES = $(addprefix ${QEDIR}/PW/src/,${PW_SOURCES})

# QE_FILES is the final list to wrap
QE_FILES = ${MODULES_FILES} ${PW_FILES}
