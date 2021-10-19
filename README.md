# QEpy - Quantum ESPRESSO in Python
   `QEpy` turns Quantum ESPRESSO (QE) into a Python DFT engine for nonstandard workflows. 

   Check out a [YouTube video](https://www.youtube.com/watch?v=cWt0BVQs-_U) with additional information (installation and examples).
   
# Contributors and funding

 - [The Quantum-Multiscale collaboration](http://www.quantum-multiscale.org/)
 - Main author: [Xuecheng Shao](mailto:xuecheng.shao@rutgers.edu) (Rutgers) 
 - Oliviero Andreussi (UNT), Davide Ceresoli (CNR, Italy), Matthew Truscott (UNT), Andrew Baczewski (Sandia), Quinn Campbell (Sandia), Michele Pavanello (Rutgers)

# Thanks to ...
 - The Quantum ESPRESSO developers for the QE codebase
 - NSF for funding the Quantum-Multiscale collaboration

# What's in this code?
Small modifications to QE routines and a quick compilation with Python wrappers. `QEpy` is based on QE-6.5 and is kept up to date with the latest QE stable release.

## Requirements
 - [Python](https://www.python.org/) (>=3.6)
 - [NumPy](https://docs.scipy.org/doc/numpy/reference/) (>=1.18.0)
 - [f90wrap](https://github.com/jameskermode/f90wrap) (latest)
 - [Quantum ESPRESSO ](https://gitlab.com/QEF/q-e/-/releases/qe-6.5) (=6.5)
 - Compiler ([GNU](https://gcc.gnu.org/fortran/)(Recommended) or [Intel](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/fortran-compiler.html))

## Optional (highly recommended):
 - [mpi4py](https://bitbucket.org/mpi4py/mpi4py) MPI for python

## Install
 - **QE**

	All static libraries should be compiled with the `-fPIC` compuiler option. Add `-fPIC` to the configuration options. E.g.,

     ```shell
	 ./configure F77=ifort F90=ifort CC=icc \
	   --with-scalapack=no -enable-openmp=yes -enable-parallel=no \
	   CFLAGS=-fPIC FFLAGS=-fPIC try_foxflags=-fPIC
     ```

	Parallel version:


     ```shell
	 ./configure F77=ifort F90=ifort CC=icc MPIF90=mpiifort \
	   --with-scalapack=no -enable-openmp=no -enable-parallel=yes \
	   CFLAGS=-fPIC FFLAGS=-fPIC try_foxflags=-fPIC
	 ```

   + `make pw` or `make pwlibs`.

 - **QEpy**

   + `git clone --recurse-submodules https://gitlab.com/shaoxc/qepy.git`
   + `qedir=${QE} python -m pip install ./qepy`

## Tips
 - `qedir` should be the folder of `QE`, which contains the *make.inc* file. This can be omitted only when the *qepy* is under the `${QE}`.
 - If not clone the submodules in the beginning, can update through `git submodule update --init --recursive`.
 - Set the *variables* can help you customize your build.

	e.g.

	- "`oldxml=yes`" can read old version QE xml file (i.e., qe-5.x).
	- "`original=yes`" only wrap original QE files and a ``qepy_mod``, which support most of versions QE.

 - If you struggle with original f90wrap, try our own modified version of f90wrap:

    ```shell
	pip install git+https://github.com/shaoxc/f90wrap.git
    ```

## FAQ
 - Some Intel MPI/MKL errors occur. What do I do?
	+ Try `export LD_PRELOAD=/opt/intel/mkl/lib/intel64/libmkl_rt.so`

 - What is the *oldxml*?
	+ The old format (`-D__OLDXML`) has been deprecated since [version 6.4](https://gitlab.com/QEF/q-e/-/releases/qe-6.4). *oldxml* allows you to read the output (wavefunctions, etc) from an old XML file.

 - How do I know when I need *oldxml*?
	+ For newer version QE (>6.3), the XML file name is 'data-file-schema.xml'. In contrast, the XML file name is 'data-file.xml' for older version QE (some with precompiler flag `-D__OLDXML`)

 - Why can't I read the wavefunctions? Why does it hang?
	+ If the wavefunctions were stored by old version QE (<6.4) or [eQE](http://eqe.rutgers.edu), please try *oldxml*.
	+ There are two different ways to store wavefunctions in QE, which is controls by PW parameter [`wf_collect`](http://www.quantum-espresso.org/Doc/INPUT_PW.html#idm68). In doubt, simply use one processor to read.

## Todo
 - ~~Update the Makefile to support Gfortran compiler.~~
 - ~~Support the updating the cell lattice.~~
 - Support hybrid functionals in embedding/iterative way.

## Bugs
 - Intel Compiler

	If you met any problems like the following, please try latest Intel compiler or GNU compiler.

	+ *[MPID_nem_tmi_pending_ssend_dequeue]: ERROR: can not find matching ssend...*
	+ The initial density totally wrong with more than one nodes.

 - OpenMPI

	If you met some problems like the following:

	+ *mca_base_component_repository_open: unable to open mca_patcher_overwrite...*

	Please update to latest version of OpenMPI, or fix with `patchelf` ([openmpi=2.1.1](https://github.com/open-mpi/ompi/issues/3705)):

    ```shell
	#!/bin/sh
	prefix="/usr/lib/x86_64-linux-gnu/openmpi"
	for filename in $(ls $prefix/lib/openmpi/*.so); do
		patchelf --add-needed libmpi.so.20 $filename
		patchelf --set-rpath "\$ORIGIN/.." $filename
	done
    ```
