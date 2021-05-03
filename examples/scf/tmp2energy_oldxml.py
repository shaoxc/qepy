import qepy
try:
    from mpi4py import MPI
except Exception:
    comm = None
else:
    comm = MPI.COMM_WORLD
    comm = comm.py2f()

fname = 'qe_in.in'
qepy.qepy_pwscf(fname, comm, True)

embed = qepy.qepy_common.embed_base()

qepy.oldxml_pw_restart.pw_readfile('header')
qepy.oldxml_pw_restart.pw_readfile('reset')
qepy.oldxml_pw_restart.pw_readfile('dim')
qepy.oldxml_pw_restart.pw_readfile('bs')
# qepy.oldxml_potinit(starting = 'file')
# qepy.oldxml_wfcinit(starting = 'file')

energy = qepy.qepy_calc_energies(embed)
qepy.qepy_stop_run(0, what = 'no')
