# QEPY - Quantum ESPRESSO Python interface
   `QEPy` turns Quantum ESPRESSO (QE) into a DFT engine for embedding or for any other purpose. 
   
# What's in this code?
Small modifications to QE routines and a quick compilation with Python wrappers. `QEPy` is based on QE-6.5 and is kept up to date with the latest QE stable release.

## Requirements
 - [Python](https://www.python.org/) (>=3.6)
 - [NumPy](https://docs.scipy.org/doc/numpy/reference/) (>=1.18.0)
 - [f90wrap](https://github.com/jameskermode/f90wrap) (latest)
 - [Quantum ESPRESSO ](https://gitlab.com/QEF/q-e/-/releases/qe-6.5) (=6.5)
 - [Intel Compiler](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/fortran-compiler.html)

## Optional (highly recommended):
 - [mpi4py](https://bitbucket.org/mpi4py/mpi4py) MPI for python

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

   + `qedir=${QE} python setup.py install --user`

## Tips
 - `qedir` should be the folder of `QE`, which contains the *make.inc* file. This can be omitted only when the *qepy* is under the `${QE}`.
 - If already build once, please delete the *build* or simply through `python setup.py clean --all`.
 - Set the *variables* can help you customize your build.

	e.g.

	- "`oldxml=yes`" can read old version QE xml file (i.e., qe-5.x).

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
 - Update the Makefile to support Gfortran compiler
 - Write a python script that can automatically update the *qepy* code according to the new version of the *QE*
