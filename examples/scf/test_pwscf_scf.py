import numpy as np
import qepy

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    comm = comm.py2f()
except Exception:
    comm = None

fname = 'qe_in.in'

qepy.qepy_pwscf(fname, comm)

embed = qepy.qepy_common.embed_base()
# embed.ldescf = True # add scf correction energy
embed.extene = 0.0
#
embed.iterative = True

for i in range(60):
    embed.mix_coef = -1.0
    qepy.qepy_electrons_scf(2, 0, embed)
    embed.mix_coef = 0.7
    qepy.qepy_electrons_scf(2, 0, embed)
    if qepy.control_flags.get_conv_elec() : break

# qepy.qepy_calc_energies(embed)
etotal = embed.etotal

#
nat = qepy.ions_base.get_nat()
forces = np.zeros((3, nat), order = 'F')
qepy.qepy_mod.qepy_set_extforces(embed, forces)
#
qepy.qepy_forces(embed = embed)
forces = qepy.force_mod.get_array_force().T

#
embed.extstress = 0.0
#
stress = np.ones((3, 3), order='F')
qepy.qepy_stress(stress, embed=embed)

qepy.punch('all')
qepy.qepy_stop_run(0, what = 'no')
