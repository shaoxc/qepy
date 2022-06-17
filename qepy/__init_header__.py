#-----------------------------------------------------------------------
#Fixed the MPI_IN_PLACE
import sys
if 'mpi4py' in sys.modules :
    import os
    import mpi4py
    from ctypes.util import find_library
    try:
        mpi4py.profile(find_library('mpi'), path = os.environ.get('LD_LIBRARY_PATH', '').split(':'))
    except Exception :
        pass
#-----------------------------------------------------------------------
