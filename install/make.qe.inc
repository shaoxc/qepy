QEDIR = $(or ${qedir}, ../../../)

include ${QEDIR}/make.inc

FOLDERS := Modules PW/src PP/src FFTXlib/src upflib LAXlib XClib
QE_FOLDERS := $(addprefix ${QEDIR}/,${FOLDERS})

QE_SOURCES += constants.f90 cell_base.f90 ions_base.f90 wavefunctions.f90 \
			  recvec.f90 control_flags.f90 io_global.f90 check_stop.f90 \
			  deviatoric.f90 qexsd.f90 environment.f90 read_input.f90 \
			  mp_pools.f90 mp_bands.f90 mp_global.f90 mp_world.f90 # Modules
QE_SOURCES += uspp.f90 uspp_param.f90 # upflib
QE_SOURCES += pwcom.f90 scf_mod.f90 punch.f90 \
			  atomic_wfc_mod.f90 close_files.f90 stress.f90 electrons.f90 \
			  scale_h.f90 pw2casino.f90 forces.f90 update_pot.f90 move_ions.f90 \
			  add_qexsd_step.f90 hinit1.f90 run_pwscf.f90 stop_run.f90 \
			  sum_band.f90 non_scf.f90 v_of_rho.f90 # PW
QE_SOURCES += mp_diag.f90 la_helper.f90 # LAXlib
QE_SOURCES += elf.f90 local_dos.f90 local_dos_mag.f90 # PP
QE_SOURCES += fft_types.f90 # FFTXlib
QE_SOURCES += dft_setting_routines.f90 # XClib

# QE_FILES is the final list to wrap
QE_FILES := $(foreach d, ${QE_FOLDERS}, $(foreach f, ${QE_SOURCES}, $(wildcard ${d}/${f})))
DIFF := $(filter-out $(notdir $(QE_FILES)), $(QE_SOURCES))
ifneq ($(DIFF),)
$(error ${DIFF} not found)
endif
