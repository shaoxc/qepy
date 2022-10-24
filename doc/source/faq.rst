.. _faq:


==========================
Frequently Asked Questions
==========================

Installation
============

What is the `oldxml`?
---------------------

  The old format (**-D__OLDXML**) has been deprecated since `version 6.4 <https://gitlab.com/QEF/q-e/-/releases/qe-6.4>`__. ``oldxml`` allows you to read the output (wavefunctions, etc) from an old XML file.

How do I know when I need `oldxml`?
-----------------------------------

  For newer version QE (>6.3), the XML file name is **data-file-schema.xml**. In contrast, the XML file name is **data-file.xml** for older version QE (some with precompiler flag **-D__OLDXML**)

Running
=======

Some Intel MPI/MKL errors occur. What do I do?
----------------------------------------------

  Try:

   .. code:: shell

    export LD_PRELOAD=/opt/intel/mkl/lib/intel64/libmkl_rt.so

Why can’t I read the wavefunctions? Why does it hang?
-----------------------------------------------------

-  If the wavefunctions were stored by old version QE (<6.4) or `eQE <http://eqe.rutgers.edu>`__, please try reinstall the QEpy with ``oldxml=yes``.

-  There are two different ways to store wavefunctions in QE, which is controls by PW parameter `wf_collect <http://www.quantum-espresso.org/Doc/INPUT_PW.html#idm68>`__.  In doubt, simply use one processor to read.

Some f90wrap errors occur. What do I do?
----------------------------------------

  + *No module named f90wrap.__main__...*

  Please try install the latest `f90wrap <https://github.com/jameskermode/f90wrap>`__ . e.g.

   .. code:: shell

    python -m pip install --upgrade f90wrap

Signal: Segmentation fault...
-----------------------------

-  `QEpy` and `mpi4py` should installed with same compiler. If not, try reinstall the `mpi4py` by:

   .. code:: shell

    python -m pip install mpi4py --no-cache-dir


Bugs
====

GCC
---
   For new version of gcc (eg. gcc=11.2) with LTO supported or under conda environment, if you met the problem when you install the QEpy_:

       +  *lto1: fatal error:...*

   Try disable link-time-optimization by add ``-fno-lto`` to the variable **FFLAGS** and **CFLAGS**.


Gfortran
--------

   #. For gfortran>=10.0 sometimes still works, but will has error:

      -  *Type mismatch between actual argument...*

      Try to add ``-fallow-argument-mismatch`` to the variable **FFLAGS** (e.g. ``FFLAGS='-fPIC -fallow-argument-mismatch'``).

   #. For some versions of the MacOS, maybe you will has error:

      - *Illegal Instruction: 4...*

      Try to add ``-mmacosx-version-min=10.14`` to the **FFLAGS**. 

     

   #. For a few versions of BLAS library, will raise error:

      -  *Segmentation fault - invalid memory reference...*

      This is due to the `zdotc` function of external libraries. More details to `here <https://gitlab.com/QEF/q-e/-/wikis/Support/zdotc-crash>`__. One solution is manually append the ``-Dzdotc=zdotc_wrapper`` to **DFLAGS** or **MANUAL_DFLAGS** in **make.inc** of QE. You also can do it during make:

     .. code:: shell

        make pw MANUAL_DFLAGS='-Dzdotc=zdotc_wrapper'


Intel Compiler
--------------

   #. If you met any problems like the following, please try a newer Intel compiler or GNU compiler.

       +  *[MPID_nem_tmi_pending_ssend_dequeue]: ERROR: can not find matching ssend...*
       +  The initial density totally wrong with more than one nodes.

   #. The gcc version between 4.8-9.2 are supported by intel compiler, which upgraded until 2022.1 version. More details to `here <https://community.intel.com/t5/Intel-oneAPI-Data-Parallel-C/Compilation-issues-with-ICPC-2021-4-and-C-14/td-p/1318571>`__.

      + *...error: attribute "__malloc__" does not take arguments...*

OpenMPI
-------

   If you met some problems like the following:

   -  *mca_base_component_repository_open: unable to open
      mca_patcher_overwrite...*

   Please update to latest version of OpenMPI, or fix with ``patchelf``
   (`openmpi=2.1.1 <https://github.com/open-mpi/ompi/issues/3705>`__):

   .. code:: shell

      #!/bin/sh
      prefix="/usr/lib/x86_64-linux-gnu/openmpi"
      for filename in $(ls $prefix/lib/openmpi/*.so); do
          patchelf --add-needed libmpi.so.20 $filename
          patchelf --set-rpath "\$ORIGIN/.." $filename
      done


.. _QEpy: https://gitlab.com/shaoxc/qepy
.. _DFTpy: http://dftpy.rutgers.edu
