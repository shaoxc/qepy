from mpi4py import MPI
import numpy as np

def comm2comm(comm):
    if comm is None:
        result=comm
    elif hasattr(comm, 'py2f'):
        result = comm.py2f()
    else:
        result = MPI.Comm.f2py(comm)
    return result

def get_data_type(data, comm=None, **kwargs):
    if comm is None: comm = MPI.COMM_WORLD
    if data is not None:
        shape = np.asarray(data.shape)
        dtype = str(data.dtype)
    else:
        shape = None
        dtype = ''
    shape = comm.bcast(shape)
    dtype = comm.bcast(dtype)
    return shape, dtype
