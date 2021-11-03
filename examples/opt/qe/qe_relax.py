#!/usr/bin/env python
# coding: utf-8

import numpy as np
import qepy

try:
    from mpi4py import MPI
    comm=MPI.COMM_WORLD.py2f()
except Exception :
    comm=None

# qepy.qepy_mod.qepy_set_stdout('qepy.out')

qepy.qepy_pwscf('vcrelax.in',comm)
embed = qepy.qepy_common.embed_base()

nstep = qepy.control_flags.get_nstep()
#-------------------------------------------------------------------------------
lforce = qepy.force_mod.get_lforce()
lstres = qepy.force_mod.get_lstres()
sigma = qepy.force_mod.get_array_sigma()
lmovecell = qepy.cellmd.get_lmovecell()

fix_volume = qepy.cell_base.get_fix_volume()
fix_area = qepy.cell_base.get_fix_area()
treinit_gvecs = qepy.control_flags.get_treinit_gvecs()

ions_status = np.ones(1, dtype = 'int32')*3
#-------------------------------------------------------------------------------
for idone in range(0, nstep):
    if idone>0 : qepy.control_flags.set_ethr(1E-6)
    qepy.electrons()

    conv_ions = True

    if lforce : qepy.forces()
    if lstres : qepy.stress(sigma)

    lmd = qepy.control_flags.get_lmd()
    lbfgs = qepy.control_flags.get_lbfgs()

    if lmd or lbfgs :
        if fix_volume: qepy.impose_deviatoric_stress(sigma)
        if fix_area: qepy.impose_deviatoric_stress_2d(sigma)

        qepy.extrapolation.update_file()

        qepy.move_ions(idone, ions_status)

        conv_ions = (ions_status == 0) or (ions_status == 1 and treinit_gvecs)

        if idone < nstep and not conv_ions :
            qepy.punch('config-nowf')

        if qepy.funct.dft_is_hybrid() : qepy.funct.stop_exx()

    if conv_ions: break

    if lmd or lbfgs :
        if ions_status == 1 :
            qepy.control_flags.set_lbfgs(False)
            qepy.control_flags.set_lmd(False)
            print('Final scf calculation at the relaxed structure', flush=True)
            qepy.reset_gvectors()
        elif ions_status == 2 :
            qepy.reset_magn()
        else :
            if treinit_gvecs :
                if lmovecell : qepy.scale_h()
                qepy.reset_gvectors()
            else :
                qepy.extrapolation.update_pot()
                qepy.hinit1()

qepy.qepy_stop_run(0, what = 'no')
