import qepy
try:
    from mpi4py import MPI
except Exception:
    comm = None
else:
    comm = MPI.COMM_WORLD
    comm = comm.py2f()

inputobj = qepy.qepy_common.input_base()
#-----------------------------------------------------------------------
inputobj.prefix = 'al'
inputobj.tmp_dir = './al.wfx/'
#-----------------------------------------------------------------------
if comm : inputobj.my_world_comm = comm
qepy.qepy_initial(inputobj)
qepy.read_file()

embed = qepy.qepy_common.embed_base()
energy = qepy.qepy_calc_energies(embed)
qepy.qepy_stop_run(0)
