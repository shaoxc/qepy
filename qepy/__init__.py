import sys
from pathlib import Path
from importlib import import_module
from importlib.metadata import version # python >= 3.8
from .__config__ import CONFIG, show as show_config
from .core import QEpyMods, find_library, load_library

__all__ = ["show_config"]

# fix MPI_IN_PLACE
# mpilib = find_library('mpifort') or find_library('mpi')
# if mpilib : load_library(mpilib)
# fix mkl only if mkl library was linked
if 'mkl_' in CONFIG.get('QE', {}).get('link', ''): load_library('mkl_rt')

path = Path(__file__).resolve().parent/'qepylibs'
sys.path.insert(0, str(path))
__path__.append(str(path))

__version__ = "7.2.0rc2"

try:
    __version__ = version("qepy")
except Exception :
    pass

# compatible with older versions API
constants = QEpyMods('qepy_modules.constants')
qepy_common = QEpyMods('qepy_pw.qepy_common')

# TODO: replace the package by subpackage
def __getattr__(pname):
    if pname.startswith('qepy_'):
        p = import_module(pname)
        return p
