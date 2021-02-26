# import atexit

# def pwscf_finalise():
#     pwpy_pwscf_finalise()


# atexit.register(pwscf_finalise)

import pkgutil
import operator
def pwpy_clean_saved():
    mods = [name for _, name, _ in pkgutil.iter_modules(pwscfpy.__path__)]
    for mod in mods :
        if hasattr(pwscfpy, mod):
            mod += '._arrays'
            operator.attrgetter(mod)(pwscfpy).clear()


pwpy_clean_saved()
