import os
import types
import sys
from functools import wraps
import pkgutil
import operator
from importlib import import_module

class Logger(type):
    def __new__(cls, name, bases, attrs):

        for attr_name, attr_value in attrs.items():
            if isinstance(attr_value, types.FunctionType):
                attrs[attr_name] = cls.stdout2file(attr_value)

        return super().__new__(cls, name, bases, attrs)

    @staticmethod
    def stdout2file(function, fileobj = None):
        return stdout2file(function, fileobj=fileobj)

    @staticmethod
    def print2file(fileobj = None):
        return print2file(fileobj=fileobj)

def print2file(fileobj = None):
    def decorator(function):
        return stdout2file(function, fileobj=fileobj)
    return decorator

def stdout2file(function, fileobj = None):
    @wraps(function)
    def wrapper(*args, **kwargs):
        fobj = fileobj
        if not fobj :
            if len(args) > 0 and hasattr(args[0], 'fileobj'): fobj = args[0].fileobj
        stdout = None
        if fobj is not None and not fobj.closed :
            if os.fstat(1).st_ino != os.fstat(fobj.fileno()).st_ino :
                stdout = os.dup(1)
                os.dup2(fobj.fileno(), 1)
        results = function(*args, **kwargs)
        if stdout is not None :
            os.dup2(stdout, 1)
            os.close(stdout)
        return results
    return wrapper

def qepy_clean_saved(module=None):
    if module is None:
        l = [v for key, v in sys.modules.items() if key.startswith('qepy_') and '.' not in key]
        for module in l:
            qepy_clean_saved(module=module)
    mods = [name for _, name, _ in pkgutil.iter_modules(module.__path__)]
    for mod in mods :
        if hasattr(module, mod):
            for item in ['_arrays', '_objs'] :
                if hasattr(operator.attrgetter(mod)(module), item):
                    attr = mod + '.' + item
                    operator.attrgetter(attr)(module).clear()


env = {
    'LOGFILE' : None, # The screen output of QE.
    'STDOUT' : None,  # file descriptor of output
    'STDIN' : None,  # file descriptor of input
    'STDIN_SAVE' : None,  # The saved screen input of QE.
    'DRIVER' : None, # save the instance of driver class
}

qepylibs = [
    'qepy_pw',
    'qepy_atomic',
    'qepy_gww_gww',
    'qepy_gww_simple_ip',
    'qepy_lr_modules',
    'qepy_pp',
    'qepy_xclib',
    'qepy_cetddft',
    'qepy_gww_head',
    'qepy_hp',
    'qepy_mbd',
    'qepy_xspectra',
    'qepy_cpv',
    'qepy_gww_minpack',
    'qepy_kcw',
    'qepy_modules',
    'qepy_pwcond',
    'qepy_dft_d3',
    'qepy_gww_pw4gww',
    'qepy_kcw_pp',
    'qepy_neb',
    'qepy_tddfpt',
    'qepy_fftxlib',
    'qepy_gww_simple',
    'qepy_ks_solvers',
    'qepy_phonon_gamma',
    'qepy_upflib',
    'qepy_gww_bse',
    'qepy_gww_simple_bse',
    'qepy_laxlib',
    'qepy_phonon_ph',
    'qepy_utilxlib',
]

class QEpyLibs(type):
    def __getattr__(cls, attr):
        if attr.startswith('qepy_'):
            # TODO: replace the package by subpackage
            p = import_module(attr)
            return p
        else :
            return object.__getattribute__(cls, attr)

class QEpyMods:
    def __init__(self, qepymod):
        self.qepymod = qepymod

    def __getattr__(self, attr):
        if attr == 'qepymod' :
            return object.__getattribute__(self, attr)
        else :
            for _, mod in self.qepymod.items() :
                if hasattr(mod, attr):
                    return getattr(mod, attr)
            else:
                raise AttributeError(f"There is no '{attr}' in {list(self.qepymod.keys())}")

    @property
    def qepymod(self):
        if isinstance(self._qepymod, str):
            self._qepymod = {self._qepymod: None}
        elif not isinstance(self._qepymod, dict):
            self._qepymod = dict.fromkeys(self._qepymod)
        for k, v in self._qepymod.items():
            if v is None :
                mod = import_module(k)
                self._qepymod[mod.__name__] = mod
        return self._qepymod

    @qepymod.setter
    def qepymod(self, value):
        self._qepymod = value
