# import atexit

# def pwscf_finalise():
#     qepy_pwscf_finalise()


# atexit.register(pwscf_finalise)

import pkgutil
import operator
def qepy_clean_saved():
    mods = [name for _, name, _ in pkgutil.iter_modules(qepy.__path__)]
    for mod in mods :
        if hasattr(qepy, mod):
            mod += '._arrays'
            operator.attrgetter(mod)(qepy).clear()


qepy_clean_saved()
__author__ = "Pavanello Research Group"
__contact__ = "m.pavanello@rutgers.edu"
__version__ = "0.0.1"
__license__ = "GPL"
__date__ = "2021-09-30"

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("qepy")
except PackageNotFoundError:
    pass
