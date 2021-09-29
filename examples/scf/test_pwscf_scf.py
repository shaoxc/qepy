import numpy as np
import qepy

try:
    from mpi4py import MPI
except Exception:
    comm = None
else:
    comm = MPI.COMM_WORLD
    comm = comm.py2f()

fname = 'qe_scf.in'

qepy.qepy_pwscf(fname, comm)

embed = qepy.qepy_common.embed_base()
embed.exttype = 0
embed.finish = False

nstep = 50
for i in range(nstep):
    if i == 0 :
        initial = True
    else :
        initial = False
    embed.initial = initial
    embed.mix_coef = -1.0
    qepy.qepy_electrons_scf(0, 0, embed)
    embed.mix_coef = 0.7
    qepy.qepy_electrons_scf(2, 0, embed)
embed.finish = True
qepy.qepy_electrons_scf(2, 0, embed)

qepy.qepy_calc_energies(embed)
etotal = embed.etotal

qepy.qepy_forces(0)
forces = qepy.force_mod.get_array_force().T

stress = np.ones((3, 3), order='F')
qepy.stress(stress)

qepy.punch('all')
qepy.qepy_stop_run(0, what = 'no')
