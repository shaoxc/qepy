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
