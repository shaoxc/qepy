.. _download_and_install:

============
Installation
============

Requirements
============

-  `Python <https://www.python.org/>`__ (>=3.6)
-  `NumPy <https://docs.scipy.org/doc/numpy/reference/>`__ (>=1.18.0)
-  `f90wrap <https://github.com/jameskermode/f90wrap>`__ (latest)
-  `Quantum ESPRESSO <https://gitlab.com/QEF/q-e/-/releases/qe-7.2>`__
   (==7.2)
-  Compiler (`GNU <https://gcc.gnu.org/fortran/>`__\ (Recommended) or
   `Intel <https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/fortran-compiler.html>`__)

Optional (highly recommended):
------------------------------

-  `mpi4py <https://mpi4py.readthedocs.io/en/stable/index.html>`__ MPI for python
-  `ASE <https://wiki.fysik.dtu.dk/ase/index.html>`__ for atomic simulation
-  `DFTpy <http://dftpy.rutgers.edu>`__ for density functional


Installation
============

Pip
---

Using pip can easy install the release version (serial) of QEpy from `PyPI <https://pypi.org/project/qepy>`_::

    $ python -m pip install qepy

.. note::

    Install the QEpy using ``pip`` only support serial version. If want to run parallel version or use some custom features, please install ``QEpy`` from source.

Source
------
    
You can download the ``QEpy`` source file from `gitlab <https://gitlab.com/shaoxc/qepy>`__.

   .. code:: shell

      git clone --recurse-submodules https://gitlab.com/shaoxc/qepy.git
      python -m pip install -U ./qepy
   
or with all features:

   .. code:: shell

      oldxml=yes ldau=yes tddft=yes python -m pip install -U ./qepy


Example on Ubuntu 22.04
+++++++++++++++++++++++

.. code:: shell


  sudo apt-get update
  sudo apt-get install --upgrade -y make git python3-dev python3-pip wget
  sudo apt-get install --upgrade -y gcc gfortran libblas-dev liblapack-dev libopenmpi-dev libfftw3-dev
  git clone --depth=1 -b qe-7.2 https://gitlab.com/QEF/q-e.git
  cd q-e
  ./configure CFLAGS=-fPIC FFLAGS='-fPIC -fallow-argument-mismatch' try_foxflags=-fPIC MPIF90=mpif90 --with-scalapack=no BLAS_LIBS='-lblas' LAPACK_LIBS='-llapack'
  make all -j 8
  make all -j 8
  cd ..
  git clone --recurse-submodules https://gitlab.com/shaoxc/qepy.git
  qedir=`pwd`/q-e/ python3 -m pip install -U ./qepy



Tips
====

-  Environment variable ``qedir`` should be the absolute path of ``QE``, which contains the *make.inc* file.
   If not set ``qedir``, the installation will download the QE code from gitlab and automatically compile it.

-  If not clone the submodules in the beginning, can update through:

   .. code:: shell

      git submodule update --init --recursive

-  Set the *variables* can help you customize your build.

   e.g.

   -  “``tddft=yes``” support real-time TDDFT by using `ce-tddft <https://github.com/dceresoli/ce-tddft>`__.

Install the QE
--------------

   The ``QE`` should be compiled before ``QEpy`` with the ``-fPIC`` (position-independent code) compiler
   option. Add ``-fPIC`` to the configuration options. E.g.,

   .. code:: shell

      ./configure CFLAGS=-fPIC FFLAGS=-fPIC try_foxflags=-fPIC

   Then,

   .. code:: shell

      make all
      export qedir=`pwd`


.. note::

    Because ``QEpy`` still under active development, non-backward-compatible changes can happen at any time. Please, clone the lastest release often.
