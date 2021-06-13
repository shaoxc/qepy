import numpy as np
import qepy

try:
    from mpi4py import MPI
except Exception:
    comm = None
else:
    comm = MPI.COMM_WORLD
    comm = comm.py2f()

fname = 'qe_in.in'

qepy.qepy_pwscf(fname, comm)

embed = qepy.qepy_common.embed_base()

qepy.qepy_electrons_scf(0, 0, embed)

qepy.qepy_calc_energies(embed)
etotal = embed.etotal

qepy.qepy_forces(0)
forces = qepy.force_mod.get_array_force().T

qepy.punch('all')
qepy.qepy_stop_run(0, what = 'no')
