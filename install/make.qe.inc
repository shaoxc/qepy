QEDIR = $(or ${qedir}, ../../../)

include ${QEDIR}/make.inc

MODULES_SOURCES = constants.f90 cell_base.f90 ions_base.f90 wavefunctions.f90 funct.f90 \
				  recvec.f90 control_flags.f90 io_global.f90 check_stop.f90 \
				  deviatoric.f90 qexsd.f90
MODULES_FILES = $(addprefix ${QEDIR}/Modules/,${MODULES_SOURCES})

PW_SOURCES = pwcom.f90 scf_mod.f90 read_file_new.f90 punch.f90 \
			 atomic_wfc_mod.f90 close_files.f90 stress.f90 electrons.f90 \
			 scale_h.f90 pw2casino.f90 forces.f90 update_pot.f90 move_ions.f90 \
			 add_qexsd_step.f90 hinit1.f90 run_pwscf.f90

PW_FILES = $(addprefix ${QEDIR}/PW/src/,${PW_SOURCES})

# QE_FILES is the final list to wrap
QE_FILES = ${MODULES_FILES} ${PW_FILES}
