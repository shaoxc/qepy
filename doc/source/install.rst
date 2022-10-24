.. _download_and_install:

============
Installation
============

Requirements
============

-  `Python <https://www.python.org/>`__ (>=3.6)
-  `NumPy <https://docs.scipy.org/doc/numpy/reference/>`__ (>=1.18.0)
-  `f90wrap <https://github.com/jameskermode/f90wrap>`__ (latest)
-  `Quantum ESPRESSO <https://gitlab.com/QEF/q-e/-/releases/qe-6.5>`__
   (=6.5)
-  Compiler (`GNU <https://gcc.gnu.org/fortran/>`__\ (Recommended) or
   `Intel <https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/fortran-compiler.html>`__)

Optional (highly recommended):
------------------------------

-  `mpi4py <https://bitbucket.org/mpi4py/mpi4py>`__ MPI for python


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


Example on Ubuntu 20.04
+++++++++++++++++++++++

.. code:: shell

  sudo apt-get update
  sudo apt-get install --upgrade make git python3-dev python3-pip
  sudo apt-get install --upgrade gcc gfortran libblas-dev liblapack-dev libopenmpi-dev libfftw3-dev
  wget https://gitlab.com/QEF/q-e/-/archive/qe-6.5/q-e-qe-6.5.tar.gz
  tar -xzvf q-e-qe-6.5.tar.gz
  cd q-e-qe-6.5
  ./configure CFLAGS=-fPIC FFLAGS=-fPIC try_foxflags=-fPIC MPIF90=mpif90 --with-scalapack=no BLAS_LIBS='-lblas' LAPACK_LIBS='-llapack'
  make pw -j 4
  cd ..
  git clone --recurse-submodules https://gitlab.com/shaoxc/qepy.git
  qedir=`pwd`/q-e-qe-6.5/ python -m pip install -U ./qepy



Tips
====

-  Environment variable ``qedir`` should be the absolute path of ``QE``, which contains the *make.inc* file.
   If not set ``qedir``, the installation will download the QE code from gitlab and automatically compile it.

-  If not clone the submodules in the beginning, can update through:

   .. code:: shell

      git submodule update --init --recursive

-  Set the *variables* can help you customize your build.

   e.g.

   -  “``oldxml=yes``” can read old version QE xml file (i.e., qe-5.x).
   -  “``ldau=yes``” will generate LDA+U (DFT+U) files based on given
      `electron configuration <https://gitlab.com/shaoxc/qepy/-/tree/master/src/ldau/qepy_econf.ini>`__.
   -  “``original=yes``” only wrap original QE files and a ``qepy_mod``,
      which also can support other version of QE
      (e.g. `6.5 <https://gitlab.com/shaoxc/qepy/-/tree/master/examples/original/6.5>`__,
      `6.8-7.1 <https://gitlab.com/shaoxc/qepy/-/tree/master/examples/original/6.8>`__).

Install the QE
--------------

   The ``QE`` should be compiled before ``QEpy`` with the ``-fPIC`` (position-independent code) compiler
   option. Add ``-fPIC`` to the configuration options. E.g.,

   .. code:: shell

      ./configure CFLAGS=-fPIC FFLAGS=-fPIC try_foxflags=-fPIC MPIF90=mpif90

   Then,

   .. code:: shell

      make pw
      export qedir=`pwd`


.. note::

    Because ``QEpy`` still under active development, non-backward-compatible changes can happen at any time. Please, clone the lastest release often.
