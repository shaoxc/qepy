import os
import subprocess
import pathlib
import shutil

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

__author__ = "Pavanello Research Group"
__contact__ = "m.pavanello@rutgers.edu"
__version__ = "0.0.1"
__license__ = "GPL"
name = 'qepy'
description = "QEPY: Quantum ESPRESSO Python interface",
long_description = """ QEPY turns Quantum ESPRESSO (QE) into a DFT engine for embedding or for any other purpose."""


class MakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        topdir = './'
        build_args = []
        env = os.environ.copy()

        try:
            import multiprocessing as mp
            nprocs = mp.cpu_count()
        except ImportError:
            nprocs = 1
        build_args += ['-j', str(nprocs)]

        # if env.get('enable_parallel', 'no').lower() == 'yes' :
        #     build_args.append('mpi')

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        makefiles = list(pathlib.Path(topdir + '/Python/').glob('*')) + \
                list(pathlib.Path(topdir + '/install/').glob('*')) + \
                list(pathlib.Path(topdir + '/src/').glob('**/*.f90'))

        for f in makefiles :
            shutil.copy(f, self.build_temp)

        subprocess.check_call(['make', '-f', 'Makefile'] + build_args, cwd=self.build_temp, env = env)

        if not os.path.exists(self.build_lib):
            os.makedirs(self.build_lib)
        elif os.path.exists(self.build_lib + os.sep + name):
            shutil.rmtree(self.build_lib+ os.sep + name)
        shutil.copytree(self.build_temp + os.sep + name, self.build_lib + os.sep + name)
        for f in pathlib.Path(self.build_temp).glob('*.so'):
            shutil.copy(f, self.build_lib + os.sep)


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
            author=__author__,
            author_email=__contact__,
            license=__license__,
            long_description=long_description,
            python_requires = '>=3.6',
            install_requires=['numpy>=1.18.0', 'f90wrap>=0.2.3'],
            extras_require={
                'mpi': [
                    'mpi4py>=3.0.0',
                    ],
                },
            packages=find_packages('./'),
            ext_modules=ext_modules,
            cmdclass = {"build_ext" : MakeBuild},
            classifiers=[
                'Development Status :: 1 - Beta',
                'Intended Audience :: Science/Research',
                'Programming Language :: Python :: 3',
                'Topic :: Scientific/Engineering :: Chemistry',
                'Topic :: Scientific/Engineering :: Physics'
                ],
            )
