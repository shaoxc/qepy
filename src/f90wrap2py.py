from sysconfig import get_config_vars
from distutils.sysconfig import get_config_vars as get_config_vars2
confvars = get_config_vars()
confvars['EXT_SUFFIX'] = '.so'
confvars['SO'] = confvars['EXT_SUFFIX']
confvars2 = get_config_vars2()
confvars2['EXT_SUFFIX'] = confvars['EXT_SUFFIX']
confvars2['SO'] = confvars['EXT_SUFFIX']

from numpy.distutils.fcompiler import load_all_fcompiler_classes
load_all_fcompiler_classes()
import sys
from f90wrap.scripts.f2py_f90wrap import main
from numpy.distutils.fcompiler import fcompiler_class
from numpy.distutils.fcompiler.gnu import Gnu95FCompiler

class Gnu95FCompilerFix(Gnu95FCompiler):
    def get_flags_linker_so(self, *args, **kwargs):
        flags = super().get_flags_linker_so(*args, **kwargs)
        if '-bundle' in flags:
            flags[flags.index('-bundle')] = '-dynamiclib'
        return flags


fcompiler_class['gnu95'] = (fcompiler_class['gnu95'][0], Gnu95FCompilerFix, fcompiler_class['gnu95'][2])
sys.exit(main())
