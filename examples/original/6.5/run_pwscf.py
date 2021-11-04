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

ndiag_ = 0
world_comm = qepy.mp_world.get_world_comm()
intra_pool_comm = qepy.mp_pools.get_intra_pool_comm()
intra_bgrp_comm = qepy.mp_bands.get_intra_bgrp_comm()
inter_bgrp_comm = qepy.mp_bands.get_inter_bgrp_comm()

qepy.mp_diag.mp_start_diag(ndiag_, world_comm, intra_pool_comm, do_distr_diag_inside_bgrp_ = True)
qepy.set_mpi_comm_4_solvers(intra_pool_comm, intra_bgrp_comm, inter_bgrp_comm)
qepy.environment.environment_start('PWSCF')
qepy.read_input.read_input_file('PW', fname)

exit_status = qepy.run_pwscf()

qepy.laxlib_free_ortho_group()
qepy.stop_run(exit_status)
qepy.do_stop(exit_status)
