import numpy as np
import qepy

# from mpi4py import MPI

# comm = MPI.COMM_WORLD
# comm = comm.py2f()
comm = None

fname = 'qe_in.in'

qepy.qepy_pwscf(fname, comm)

embed = qepy.qepy_common.embed_base()

nr = np.zeros(3, dtype = 'int32')
qepy.qepy_mod.qepy_get_grid(nr)
extpot = np.zeros((np.prod(nr), 1), order = 'F')
qepy.qepy_mod.qepy_set_extpot(embed, extpot)
embed.exttype = 0

qepy.qepy_electrons_scf(0, 0, embed)

qepy.qepy_calc_energies(embed)
etotal = embed.etotal

qepy.qepy_forces(0)
forces = qepy.force_mod.get_array_force().T

rho = np.zeros((np.prod(nr), 1), order = 'F')
qepy.qepy_mod.qepy_get_rho(rho)

qepy.punch('all')
qepy.qepy_stop_run(0, what = 'no')

# print('Forces :\n', forces)
# print('Total energy :', etotal)
# from dftpy.constants import ENERGY_CONV, LEN_CONV
# from dftpy.formats import ase_io
# ions = ase_io.ase_read(fname, format = 'espresso-in')
