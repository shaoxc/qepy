import qepy

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    commf = comm.py2f()
except Exception:
    comm = None
    commf = None

fname = 'qe_in.in'
qepy.mp_global.mp_startup(my_world_comm=commf, start_images=True)
qepy.environment.environment_start('PWSCF')

qepy.read_input.read_input_file('PW', fname)

exit_status = qepy.run_pwscf()

qepy.laxlib_end()
qepy.stop_run(exit_status)
qepy.do_stop(exit_status)
