# QEPY - Quantum ESPRESSO Python interface
   `QEPy` turns Quantum ESPRESSO (QE) into a DFT engine for embedding or for any other purpose. 
   
# What's in this code?
Small modifications to QE routines and a quick compilation with Python wrappers. `QEPy` is based on QE-6.5 and is kept up to date with the latest QE stable release.

## Requirements
 - Python 3.6 or later
 - Numpy >= 1.18
 - Intel compiler (ifort)
 - f90wrap package (latest)
 - Quantum ESPRESSO source distribution (qe-6.5)

## Install
 - **f90wrap**

    Install the origin version:

    ```shell
	pip install git+https://github.com/jameskermode/f90wrap
    ```
	If you struggle, try our own modified version of f90wrap:

    ```shell
	pip install git+https://github.com/shaoxc/f90wrap
    ```



 - **QE**

	All static libraries must be compiled with the `-fPIC` compuiler option. Add `-fPIC` to the configuration options. E.g.,

     ```shell
	 ./configure F77=ifort F90=ifort CC=icc \
	   --with-scalapack=no -enable-openmp=yes -enable-parallel=no \
	   CFLAGS=-fPIC FFLAGS='-mcmodel=large -fPIC' FOXFLAGS=-fPIC
     ```

	Parallel version:


     ```shell
	 ./configure F77=ifort F90=ifort CC=icc MPIF90=mpiifort \
	   --with-scalapack=intel -enable-openmp=no -enable-parallel=yes \
	   CFLAGS=-fPIC FFLAGS='-mcmodel=large -fPIC' FOXFLAGS=-fPIC
	 ```

   + After configuration, you also need add `-fPIC` to `FOX_FLAGS` in the *make.inc* file.
   + Build the normal ***pw*** or ***pwlibs***.

 - **QEPY**

   + Copy the *qepy* into the ${QE} directory.
   + Go to *qepy* directory and `make` (serial) or `make mpi` (parallel).
   + Go to *qepy* directory and `make install`.

## Tips
 - The ***QE*** and ***QEPY*** should be both either serial or parallel.
 - `make help` will show the Makefile.
 - The *variables* of Makefile help you customize your build.

	e.g.

	- "`export oldxml=yes`" can read old version QE xml file (i.e., qe-5.x).
	- "`export prefix=~/.local/lib/python3.8/site-packages/`" can set the folder for installation.

## FAQ
 - Some Intel MPI/MKL errors occur. What do I do?
	+ Make sure QE and QEPY are of same version 
	+ Try `export LD_PRELOAD=/opt/intel/mkl/lib/intel64/libmkl_rt.so`

 - What is the *oldxml*?
	+ The old format (-D\_\_OLDXML) has been deprecated since [version 6.4](https://github.com/QEF/q-e/releases/tag/qe-6.4). *oldxml* allows you to read the output (wavefunctions, etc) from an old XML file.

 - Why can't I read the wavefunctions? Why does it hang?
	+ If the wavefunctions were stored by old version QE (<6.4) or [eQE](http://eqe.rutgers.edu), please try *oldxml*.
	+ There are two different ways to store wavefunctions in QE, which is controls by PW parameter [`wf_collect`](http://www.quantum-espresso.org/Doc/INPUT_PW.html#idm68). In doubt, simply use one processor to read.

## Todo
 - Update the Makefile to support Gfortran compiler
 - Write a python script that can automatically update the *qepy* code according to the new version of the *QE*
