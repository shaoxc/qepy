QEDIR = $(or ${qedir}, ../../../)

include ${QEDIR}/make.inc

MODULES_SOURCES = constants.f90 cell_base.f90 ions_base.f90 wavefunctions.f90 \
				  recvec.f90 control_flags.f90 io_global.f90 check_stop.f90 \
				  deviatoric.f90 qexsd.f90 environment.f90 read_input.f90 \
				  mp_pools.f90 mp_bands.f90 mp_global.f90 mp_world.f90 uspp.f90
MODULES_FILES = $(addprefix ${QEDIR}/Modules/,${MODULES_SOURCES})

PW_SOURCES = pwcom.f90 scf_mod.f90 punch.f90 \
			 atomic_wfc_mod.f90 close_files.f90 stress.f90 electrons.f90 \
			 scale_h.f90 pw2casino.f90 forces.f90 update_pot.f90 move_ions.f90 \
			 add_qexsd_step.f90 hinit1.f90 run_pwscf.f90 stop_run.f90 \
			 sum_band.f90 non_scf.f90 v_of_rho.f90
PW_FILES = $(addprefix ${QEDIR}/PW/src/,${PW_SOURCES})

LAXLIB_SOURCES = mp_diag.f90 la_helper.f90
LAXLIB_FILES = $(addprefix ${QEDIR}/LAXlib/,${LAXLIB_SOURCES})

PP_SOURCES = elf.f90 local_dos.f90 local_dos_mag.f90
PP_FILES = $(addprefix ${QEDIR}/PP/src/,${PP_SOURCES})

UTILXLIB_SOURCES = set_mpi_comm_4_solvers.f90
UTILXLIB_FILES = $(addprefix ${QEDIR}/UtilXlib/,${UTILXLIB_SOURCES})

FFTXLIB_SOURCES = fft_types.f90
FFTXLIB_FILES = $(addprefix ${QEDIR}/FFTXlib/,${FFTXLIB_SOURCES})

# QE_FILES is the final list to wrap
QE_FILES = ${MODULES_FILES} ${PW_FILES} ${LAXLIB_FILES} ${UTILXLIB_FILES} ${FFTXLIB_FILES} ${PP_FILES}
