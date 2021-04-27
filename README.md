# QEPY - Quantum ESPRESSO Python interface
   In order to use Quantum ESPRESSO in embedding code, we modified some code of quantum-espresso and compile the Python wrappers. The source code was modified based on QE-6.5 and commented out unused code.

## Requirements
 - Python 3.6 or later
 - Numpy >= 1.18
 - Intel compiler (ifort)
 - f90wrap package (latest)
 - Quantum-espresso source distribution (qe-6.5)

## Install
 - Install f90wrap:

   <!--Due to the QE was written more freely than normal code. So we need modified the f90wrap a little bit. You can directly use our modified version:-->
   <!--pip install git+https://gitlab.com/shaoxc/f90wrap@modified-->
    
    ```shell
	pip install git+https://github.com/jameskermode/f90wrap
    ```

 - Install QE

   + All static libraries must be compiled with the `-fPIC` compuiler option, so you need add `-fPIC` for all configuration. e.g.

     ```shell
	 ./configure F77=ifort F90=ifort CC=icc \
	    --with-scalapack=no -enable-openmp=yes -enable-parallel=no \
	 	CFLAGS=-fPIC FFLAGS='-mcmodel=large -fPIC' 
     ```
    or parallel version:

     ```shell
	 ./configure F77=ifort MPIF90=mpiifort F90=ifort CC=icc \
	   --with-scalapack=intel FFLAGS=-mcmodel="large -fPIC" \
	   -enable-openmp=no -enable-parallel=yes CFLAGS=-fPIC FOXFLAGS=-fPIC
	 ```

   + After configuration, you also need add `-fPIC` to `FOX_FLAGS` in the *make.inc* file.
   + Build the normal ***pw*** or ***pwlibs***.

 - Install PWSCFPY

   + Copy the *qepy* into the ${QE}/PW directory.
   + Go to *qepy* directory and `make` (serial) or `make mpi` (parallel).
   + Go to *qepy* directory and `make install`.

## TODO
 - Update the Makefile to support Gfortran compiler
 - Write a python script that can automatically update the *qepy* code according to the new version of the *QE*
