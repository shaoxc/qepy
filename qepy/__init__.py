import sys
import os
from pathlib import Path
from ctypes import util, CDLL, RTLD_LOCAL, RTLD_GLOBAL
path = Path(__file__).resolve().parent/'qepylibs'
sys.path.insert(0, str(path))
# FIX MPI_IN_PLACE and MKL
if 'mpi4py' in sys.modules :
    if hasattr(util, '_findLib_ld') and hasattr(util, '_get_soname') :
        mpilib = util._get_soname(util._findLib_ld('mpi'))
    else :
        mpilib = None
    mpilib = mpilib or util.find_library('mpi') or util.find_library('mpifort')
    try:
        CDLL(mpilib, RTLD_LOCAL | RTLD_GLOBAL)
    except Exception :
        pass
try:
    if hasattr(util, '_findLib_ld'):
        mkllib = os.path.basename(util._findLib_ld('mkl_rt'))
    else :
        mkllib = util.find_library('mkl_rt')
    CDLL(mkllib, RTLD_LOCAL | RTLD_GLOBAL)
except Exception :
    pass
# End Fix

from qepy.driver import Driver

__author__ = "Pavanello Research Group"
__contact__ = "m.pavanello@rutgers.edu"
__version__ = "7.2.0rc0"
__license__ = "GPL"
__date__ = "2023-11-30"

try:
    from importlib.metadata import version # python >= 3.8
except Exception :
    try:
        from importlib_metadata import version
    except Exception :
        pass

try:
    __version__ = version("qepy")
except Exception :
    pass
