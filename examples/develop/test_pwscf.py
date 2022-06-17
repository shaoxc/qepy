import numpy as np
import qepy

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    commf = comm.py2f()
except Exception:
    comm = None
    commf = None

fname = 'qe_in.in'

qepy.qepy_pwscf(fname, commf)

embed = qepy.qepy_common.embed_base()

qepy.qepy_electrons_scf(2, 0, embed)

nscf = qepy.control_flags.get_n_scf_steps()
conv_flag = bool(qepy.control_flags.get_conv_elec())
info = 'Converged {} at {} steps'.format(conv_flag, nscf)
qepy.qepy_mod.qepy_write_stdout(info)

# get density
nr = qepy.qepy_mod.qepy_get_grid()
nspin = qepy.lsda_mod.get_nspin()
rho = np.empty((np.prod(nr), nspin), order = 'F')
qepy.qepy_mod.qepy_get_rho(rho)
ncharge = rho.sum()*qepy.cell_base.get_omega()/np.prod(nr)
if comm is None or comm.rank == 0 : print('ncharge', ncharge, flush = True)

# qepy.qepy_calc_energies(embed)
etotal = embed.etotal

qepy.qepy_forces()
forces = qepy.force_mod.get_array_force().T

stress = np.ones((3, 3), order='F')
qepy.qepy_stress(stress)

qepy.punch('all')
qepy.qepy_stop_run(0, what = 'no')
