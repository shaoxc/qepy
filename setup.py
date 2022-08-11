import os
import sys
import re
import subprocess
import pathlib
import shutil

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

with open('qepy/__init__.py') as fd :
    lines = fd.read()
    __version__ = re.search('__version__ = "(.*)"', lines).group(1)
    __author__ = re.search('__author__ = "(.*)"', lines).group(1)
    __contact__ = re.search('__contact__ = "(.*)"', lines).group(1)
    __license__ = re.search('__license__ = "(.*)"', lines).group(1)

name = 'qepy'
description = "QEPY: Quantum ESPRESSO Python interface",
long_description = """ QEPY turns Quantum ESPRESSO (QE) into a DFT engine for embedding or for any other purpose."""

fix_mpi4py = """#Fix the MPI_IN_PLACE
import sys
if 'mpi4py' in sys.modules :
    import os
    import mpi4py
    from ctypes.util import find_library
    try:
        mpi4py.profile(find_library('mpi'), path = os.environ.get('LD_LIBRARY_PATH', '').split(':'))
    except Exception :
        pass
#End fix
"""


class MakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        topdir = './'
        build_args = []
        env = os.environ.copy()
        self.build_name = self.build_lib+ os.sep + name

        try:
            import multiprocessing as mp
            nprocs = max(mp.cpu_count()//2, 2)
        except ImportError:
            nprocs = 4
        build_args += ['-j', str(nprocs)]

        if os.path.exists(self.build_temp): shutil.rmtree(self.build_temp)
        #remove *.so files
        # for f in pathlib.Path(self.build_temp).glob('*.so'): os.remove(f)

        if not os.path.exists(self.build_temp): os.makedirs(self.build_temp)

        makefiles = list(pathlib.Path(topdir + '/qepy/').glob('__*__.py')) + \
            list(pathlib.Path(topdir + '/install/').glob('*')) + \
            list(pathlib.Path(topdir + '/src/').glob('*'))

        qeversion = 'qe-' + env.get('qeversion', '6.5')

        for f in makefiles :
            if f.is_file():
                shutil.copy2(f, self.build_temp)
            else :
                if f.name.startswith('qe-') :
                    if f.name != qeversion : continue
                shutil.copytree(f, self.build_temp+os.sep+f.name)

        makefiles = list(pathlib.Path(self.build_temp).glob('**/*.f90'))
        for f in makefiles : shutil.copy2(f, self.build_temp)

        cwd = os.getcwd()
        os.chdir(self.build_temp)
        if env.get('ldau', 'no').lower() == 'yes' :
            sys.path.insert(0, './')
            from ldau.qepy_ldau_patch import ini2files
            ini2files('ldau/qepy_econf.ini')
        os.chdir(cwd)

        env['PYTHON'] = sys.executable

        if env.get('tddft', 'no').lower() == 'yes' :
            subprocess.check_call(['make', '-f', 'Makefile.cetddft'] + build_args, cwd=self.build_temp, env = env)

        subprocess.check_call(['make', '-f', 'Makefile'] + build_args, cwd=self.build_temp, env = env)

        with open(self.build_temp + '/qepy/__init__.py', 'r+') as fh :
            lines = fh.readlines()
            fh.seek(0)
            first = True
            for line in lines :
                if first and '_qepy' in line :
                    fh.write(fix_mpi4py)
                    first = False
                fh.write(line)

        if not os.path.exists(self.build_lib): os.makedirs(self.build_lib)
        # if os.path.exists(self.build_name): shutil.rmtree(self.build_name)
        for f in pathlib.Path(self.build_temp + os.sep + name).glob('*'):
            shutil.copy2(f, self.build_name)
        for f in pathlib.Path(self.build_temp).glob('*.so'):
            shutil.copy2(f, self.build_lib + os.sep)


extensions_qepy = Extension(
        "qepy._qepy",
        sources = [],
        )

ext_modules = [extensions_qepy, ]

if __name__ == "__main__":
    setup(
            name=name,
            url='https://gitlab.com/shaoxc/qepy',
            description=description,
            version=__version__,
            use_scm_version={'version_scheme': 'post-release'},
            author=__author__,
            author_email=__contact__,
            license=__license__,
            long_description=long_description,
            python_requires = '>=3.6',
            setup_requires=[
                'setuptools_scm',
                'f90wrap @ git+https://github.com/jameskermode/f90wrap.git@master',
                ],
            install_requires=[
                'setuptools_scm',
                'numpy>=1.18.0',
                'f90wrap @ git+https://github.com/jameskermode/f90wrap.git@master',
                'importlib-metadata>=0.12;python_version<"3.8"'
                ],
            extras_require={
                'mpi': [
                    'mpi4py>=3.0.2',
                    ],
                },
            packages=find_packages('./'),
            ext_modules=ext_modules,
            # include_package_data=True,
            cmdclass = {"build_ext" : MakeBuild},
            classifiers=[
                'Development Status :: 1 - Beta',
                'Intended Audience :: Science/Research',
                'Programming Language :: Python :: 3',
                'Topic :: Scientific/Engineering :: Chemistry',
                'Topic :: Scientific/Engineering :: Physics'
                ],
            zip_safe = False,
            )
