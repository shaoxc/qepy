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

1. QE
-----

   All libraries should be compiled with the ``-fPIC`` (position-independent code) compuiler
   option. Add ``-fPIC`` to the configuration options. E.g.,

   .. code:: shell

      ./configure CFLAGS=-fPIC FFLAGS=-fPIC try_foxflags=-fPIC MPIF90=mpif90

   or intel compiler:

   .. code:: shell

      ./configure CFLAGS=-fPIC FFLAGS=-fPIC try_foxflags=-fPIC MPIF90=mpiifort

   Then,

   .. code:: shell

      make pw

2. QEpy
-------

   Installation:

   .. code:: shell

      git clone --recurse-submodules https://gitlab.com/shaoxc/qepy.git
      qedir=${QE} python -m pip install -U ./qepy
   
   or with all features:

   .. code:: shell

      qedir=${QE} oldxml=yes ldau=yes tddft=yes python -m pip install -U ./qepy



Tips
====

-  ``qedir`` should be the absolute path of ``QE``, which contains the
   *make.inc* file. This can be omitted only when the *qepy* is under
   the ``${QE}``.

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
      `6.8 <https://gitlab.com/shaoxc/qepy/-/tree/master/examples/original/6.8>`__).

-  If you struggle with original f90wrap, try our own modified version
   of f90wrap:

   .. code:: shell

      pip install git+https://github.com/shaoxc/f90wrap.git



.. note::

    Because ``QEpy`` still under active development, non-backward-compatible changes can happen at any time. Please, clone the lastest release often.
