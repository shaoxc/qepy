# PWSCFPY - PWSCF (quantum-espresso) Python interface
   In order to use PWSCF in embedding code, we modified some code of quantum-espresso and compile the Python wrappers. The source code was modified based on QE-6.5 and commented out unused code.

## Requirements
 - Python 3.6 or later
 - Numpy >= 1.18
 - Intel compiler (ifort)
 - f90wrap package (Modified)
 - Quantum-espresso source distribution (qe-6.5)

## Install
 - Install f90wrap:

   Due to the QE was written more freely than normal code. So we need modified the f90wrap a little bit. You can directly use our modified version:
    
    ```shell
	pip install git+https://gitlab.com/shaoxc/f90wrap@modified
    ```

 - Install QE

   + All static libraries must be compiled with the `-fPIC` compuiler option, so you need add `-fPIC` for all configuration. e.g.

     ```shell
	 ./configure F77=ifort MPIF90=mpiifort F90=ifort CC=icc --with-scalapack=intel \
	 				-enable-openmp=yes -enable-parallel=no \
	 				CFLAGS=-fPIC FFLAGS='-mcmodel=large -fPIC' 
     ```
   + After configuration, you also need add `-fPIC` to `FOX_FLAGS` in the *make.inc* file.
   + Build the normal ***pw*** or ***pwlibs***.

 - Install PWSCFPY

   + Copy the *pwpy* into the ${QE}/PW directory.
   + Go to *pwpy* directory and `make`.

## TODO
 - Update the Makefile to support Gfortran compiler
 - Write a python script that can automatically update the *pwpy* code according to the new version of the *QE*
