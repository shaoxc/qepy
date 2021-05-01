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

qepy.qepy_pw_restart_new.qepy_read_xml_file(alloc=False)
qepy.qepy_potinit(starting = 'file')
qepy.qepy_wfcinit(starting = 'file')

energy = qepy.qepy_calc_energies(embed)
qepy.qepy_stop_run(0, what = 'no')
