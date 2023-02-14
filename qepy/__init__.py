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
__date__ = "2023-02-14"

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
