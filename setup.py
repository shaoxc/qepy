import os
import sys
import re
import subprocess
from pathlib import Path
import shutil

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

if os.environ.get('qepybackend', 'meson').lower() != 'meson' :
    # setuptools > 59.8.0
    os.environ['SETUPTOOLS_USE_DISTUTILS'] = 'stdlib'
    if sys.version_info >= (3, 12):
        raise RuntimeError('Only meson support python>=3.12')

with open('qepy/__init__.py') as fh :
    lines = fh.read()
    __version__ = re.search('__version__ = "(.*)"', lines).group(1)
    __author__ = re.search('__author__ = "(.*)"', lines).group(1)
    __contact__ = re.search('__contact__ = "(.*)"', lines).group(1)
    __license__ = re.search('__license__ = "(.*)"', lines).group(1)

with open('README.md') as fh :
    long_description = fh.read()

name = 'qepy'
description = "QEpy: Quantum ESPRESSO Python interface"

class MakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        topdir = Path('.')
        build_args = ' '
        env = os.environ
        self.build_name = Path(self.build_lib) / name
        self.build_path = Path(self.build_temp)
        env['PYTHON'] = sys.executable

        try:
            import multiprocessing as mp
            nprocs = max(mp.cpu_count()//2, 2)
        except ImportError:
            nprocs = 4
        build_args += ' -j ' + str(nprocs)

        if env.get('qepydev', 'no').lower() == 'yes' :
            print("Keep the previous compiled files", flush = True)
        else :
            if self.build_path.is_dir(): shutil.rmtree(self.build_temp)
            shutil.copytree(topdir / 'src', self.build_temp)

        # install the QE
        qedir = env.get('qedir', '')
        if not qedir :
            qedir = self.build_temp + '/q-e'
            qe_download = "git clone -b qe-7.2 --depth=1 https://gitlab.com/QEF/q-e.git "+ qedir
            subprocess.check_call(qe_download, env = env, shell=True, text=True)
            if '-fPIC' not in env.get('CFLAGS', ''):
                env['CFLAGS'] = '-fPIC ' + env.get('CFLAGS', '')
            if '-fPIC' not in env.get('FFLAGS', ''):
                env['FFLAGS'] = '-fPIC ' + env.get('FFLAGS', '')
            if '-fPIC' not in env.get('try_foxflags', ''):
                env['try_foxflags'] = '-fPIC ' + env.get('try_foxflags', '')
            # blas and lapack
            if 'BLAS_LIBS' not in env : env['BLAS_LIBS'] = '-lblas'
            if 'LAPACK_LIBS' not in env : env['LAPACK_LIBS'] = '-llapack'
            qe_install_flags = env.get('QE_INSTALL_FLAGS', '')
            print('env', env)
            print('build_args', build_args)
            res = subprocess.run("./configure " + qe_install_flags, cwd=qedir, env = env, shell=True, capture_output=True, text=True)
            if res.returncode > 0 :
                print('Some errors happened in configure...')
                print(res.stderr)
                subprocess.run("cat install/config.log", cwd=qedir, env = env, shell=True)
                raise RuntimeError('QE configure failed.')

            res = subprocess.run("make all " + build_args, cwd=qedir, env = env, shell=True, capture_output=True, text=True)
            if res.returncode > 0 :
                print(res.stderr[-100:])
                print("'make w90' sometimes will failed at first time, so try again")
                res = subprocess.run("make all " + build_args, cwd=qedir, env = env, shell=True, capture_output=True, text=True)
                if res.returncode > 0 :
                    print(res.stderr)
                    raise RuntimeError('QE installation failed.')

            env['qedir'] = os.path.abspath(qedir)

        res = subprocess.run('make all ' + build_args, cwd=self.build_temp, env = env, shell=True, capture_output=True, text=True)
        if res.returncode > 0 :
            print("stdout:", res.stdout[-10000:])
            print("stderr:", res.stderr)
            raise RuntimeError('QEpy installation failed.')

        if env.get('tddft', 'no').lower() == 'yes' :
            res = subprocess.run('make qepy_cetddft' + build_args, cwd=self.build_temp, env = env, shell=True, capture_output=True, text=True)
            if res.returncode > 0 :
                print("stderr:", res.stderr)
                raise RuntimeError('QEpy[cetddft] installation failed.')

        qepylibs = self.build_name / 'qepylibs'
        for f in self.build_path.glob('*.so'):
            shutil.copy2(f, qepylibs)
        for f in self.build_path.glob('qepy_*'):
            if f.is_file():
                shutil.copy2(f, self.build_name)
            else :
                target = qepylibs / f.name
                if target.is_dir(): shutil.rmtree(target)
                shutil.copytree(f, target)
        #
        shutil.copy2(self.build_path / '__config__.py', self.build_name)


extensions_qepy = Extension(
        "qepy._qepy",
        sources = [],
        )

ext_modules = [extensions_qepy, ]

release = 1
if release :
    VERSION = {'version' : __version__}
else :
    VERSION = {
            'use_scm_version': {'version_scheme': 'post-release'},
            'setup_requires': ['setuptools_scm'],
            }

if __name__ == "__main__":
    setup(
            name=name,
            url='https://gitlab.com/shaoxc/qepy',
            description=description,
            author=__author__,
            author_email=__contact__,
            license=__license__,
            **VERSION,
            long_description=long_description,
            long_description_content_type='text/markdown',
            python_requires = '>=3.8',
            install_requires=[
                'setuptools',
                'numpy>=1.19.2',
                'meson>=0.63.3',
                'ninja>=1.8.2',
                'f90wrap>=0.2.14',
                # 'f90wrap@git+https://github.com/jameskermode/f90wrap',
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
                'Development Status :: 3 - Alpha',
                'Intended Audience :: Science/Research',
                'Programming Language :: Python :: 3',
                'Programming Language :: Python :: 3.8',
                'Programming Language :: Python :: 3.9',
                'Programming Language :: Python :: 3.10',
                'Programming Language :: Python :: 3.11',
                'Programming Language :: Python :: 3.12',
                'Topic :: Scientific/Engineering :: Chemistry',
                'Topic :: Scientific/Engineering :: Physics'
                ],
            zip_safe = False,
            )
