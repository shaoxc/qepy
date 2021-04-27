import pwscfpy
try:
    from mpi4py import MPI
except Exception:
    comm = None
else:
    comm = MPI.COMM_WORLD
    comm = comm.py2f()

inputobj = pwscfpy.pwpy_common.input_base()
#-----------------------------------------------------------------------
inputobj.prefix = 'al'
inputobj.tmp_dir = './al.wfx/'
#-----------------------------------------------------------------------
if comm : inputobj.my_world_comm = comm
pwscfpy.pwpy_initial(inputobj)
pwscfpy.read_file()

embed = pwscfpy.pwpy_common.embed_base()
energy = pwscfpy.pwpy_calc_energies(embed)
pwscfpy.pwpy_stop_run(0)
