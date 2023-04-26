import os
import sys
import re
import subprocess
import pathlib
import shutil

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

with open('qepy/__init__.py') as fh :
    lines = fh.read()
    __version__ = re.search('__version__ = "(.*)"', lines).group(1)
    __author__ = re.search('__author__ = "(.*)"', lines).group(1)
    __contact__ = re.search('__contact__ = "(.*)"', lines).group(1)
    __license__ = re.search('__license__ = "(.*)"', lines).group(1)

with open('qepy/__new__.py') as fh :
    init_new = fh.read()

with open('README.md') as fh :
    long_description = fh.read()

name = 'qepy'
description = "QEpy: Quantum ESPRESSO Python interface"

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

        if env.get('qepydev', False):
            print("only remove *.so files", flush = True)
            for f in pathlib.Path(self.build_temp).glob('*.so'): os.remove(f)
        else :
            if os.path.exists(self.build_temp): shutil.rmtree(self.build_temp)

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

        # install the QE
        qedir = env.get('qedir', '')
        if not qedir :
            qedir = self.build_temp + '/q-e'
            qe_download = ["git", "clone", "-b", "qe-6.5", "--depth=1", "https://gitlab.com/QEF/q-e.git", qedir]
            subprocess.check_call(qe_download, env = env)
            if '-fPIC' not in env.get('CFLAGS', ''):
                env['CFLAGS'] = '-fPIC ' + env.get('CFLAGS', '')
            if '-fPIC' not in env.get('FFLAGS', ''):
                env['FFLAGS'] = '-fPIC ' + env.get('FFLAGS', '')
            if '-fPIC' not in env.get('try_foxflags', ''):
                env['try_foxflags'] = '-fPIC ' + env.get('try_foxflags', '')
            # blas and lapack
            if 'BLAS_LIBS' not in env : env['BLAS_LIBS'] = '-lblas'
            if 'LAPACK_LIBS' not in env : env['LAPACK_LIBS'] = '-llapack'
            qe_install_flags = env.get('QE_INSTALL_FLAGS', '').split('|')
            print('env', env)
            res = subprocess.run(["./configure"] + qe_install_flags, cwd=qedir, env = env, check=False, capture_output=True)
            stderr=res.stderr.decode()
            if 'error' in stderr:
                print('Some errors happened in configure...')
                subprocess.run(["cat", "install/config.log"], cwd=qedir, env = env, check = False)
                print(stderr)
                exit()
            subprocess.check_call(["make", "pwall"] + build_args, cwd=qedir, env = env)
            env['qedir'] = os.path.abspath(qedir)

        if env.get('tddft', 'no').lower() == 'yes' :
            subprocess.check_call(['make', '-f', 'Makefile.cetddft'] + build_args, cwd=self.build_temp, env = env)

        subprocess.check_call(['make', '-f', 'Makefile'] + build_args, cwd=self.build_temp, env = env)

        with open(self.build_temp + '/qepy/__init__.py', 'r+') as fh :
            lines = fh.readlines()
            fh.seek(0)
            first = True
            for line in lines :
                if first and '_qepy' in line :
                    fh.write(init_new)
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
            # use_scm_version={'version_scheme': 'post-release'},
            author=__author__,
            author_email=__contact__,
            license=__license__,
            long_description=long_description,
            long_description_content_type='text/markdown',
            python_requires = '>=3.7',
            install_requires=[
                'setuptools_scm',
                'setuptools<=59.8.0',
                'numpy>=1.18.0',
                'f90wrap>=0.2.8',
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
                'Development Status :: 3 - Alpha',
                'Intended Audience :: Science/Research',
                'Programming Language :: Python :: 3',
                'Programming Language :: Python :: 3.7',
                'Programming Language :: Python :: 3.8',
                'Programming Language :: Python :: 3.9',
                'Programming Language :: Python :: 3.10',
                'Topic :: Scientific/Engineering :: Chemistry',
                'Topic :: Scientific/Engineering :: Physics'
                ],
            zip_safe = False,
            )
