import numpy as np
import pwscfpy

from mpi4py import MPI

comm = MPI.COMM_WORLD
comm = comm.py2f()
# comm = None

fname = 'qe_in.in'

pwscfpy.pwpy_pwscf(fname, comm)

embed = pwscfpy.pwpy_common.embed_base()

nr = np.zeros(3, dtype = 'int32')
pwscfpy.pwpy_mod.pwpy_get_grid(nr)
extpot = np.zeros((np.prod(nr), 1), order = 'F')
pwscfpy.pwpy_mod.pwpy_set_extpot(embed, extpot)
embed.exttype = 0

pwscfpy.pwpy_electrons_scf(0, 0, embed)

pwscfpy.pwpy_calc_energies(embed)
etotal = embed.etotal

pwscfpy.pwpy_forces(0)
forces = pwscfpy.force_mod.get_array_force().T

rho = np.zeros((np.prod(nr), 1), order = 'F')
pwscfpy.pwpy_mod.pwpy_get_rho(rho)
pwscfpy.pwpy_stop_run(0)

# print('Forces :\n', forces)
# print('Total energy :', etotal)
# from dftpy.constants import ENERGY_CONV, LEN_CONV
# from dftpy.formats import ase_io
# ions = ase_io.ase_read(fname, format = 'espresso-in')
