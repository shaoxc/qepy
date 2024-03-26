import sys
import os
from pathlib import Path
from ctypes import util, CDLL, RTLD_LOCAL, RTLD_GLOBAL
from importlib import import_module
from importlib.metadata import version # python >= 3.8
path = Path(__file__).resolve().parent/'qepylibs'
sys.path.insert(0, str(path))
__path__.append(str(path))
# FIX MPI_IN_PLACE and MKL
if hasattr(util, '_findLib_ld') and hasattr(util, '_get_soname') :
    mpilib = util._get_soname(util._findLib_ld('mpi'))
else :
    mpilib = None
mpilib = util.find_library('mpifort') or mpilib or util.find_library('mpi')
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

__author__ = "Pavanello Research Group"
__contact__ = "m.pavanello@rutgers.edu"
__version__ = "7.2.0rc0"
__license__ = "GPL"
__date__ = "2024-02-19"

try:
    __version__ = version("qepy")
except Exception :
    pass

# compatible with older versions API
from .core import QEpyMods
constants = QEpyMods('qepy_modules.constants')
qepy_common = QEpyMods('qepy_pw.qepy_common')

# TODO: replace the package by subpackage
def __getattr__(pname):
    if pname.startswith('qepy_'):
        p = import_module(pname)
        return p
