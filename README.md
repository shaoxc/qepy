# QEpy - Quantum ESPRESSO in Python
   `QEpy` turns Quantum ESPRESSO (QE) into a Python DFT engine for nonstandard workflows.

   <!--Check out a [YouTube video](https://www.youtube.com/watch?v=cWt0BVQs-_U) with additional information (installation and examples).-->

## Contributors and funding
 - [The Quantum-Multiscale collaboration](http://www.quantum-multiscale.org/)
 - Main author: [Xuecheng Shao](mailto:xuecheng.shao@rutgers.edu) (Rutgers)
 - Oliviero Andreussi (UNT), Davide Ceresoli (CNR, Italy), Matthew Truscott (UNT), Andrew Baczewski (Sandia), Quinn Campbell (Sandia), Michele Pavanello (Rutgers)


## Thanks to ...
 - The Quantum ESPRESSO developers for the QE codebase
 - NSF for funding the Quantum-Multiscale collaboration

## Requirements
 - [Python](https://www.python.org/) (>=3.7)
 - [NumPy](https://docs.scipy.org/doc/numpy/reference/) (>=1.18.0)
 - [f90wrap](https://github.com/jameskermode/f90wrap) (>=0.2.8)
 - [Quantum ESPRESSO ](https://gitlab.com/QEF/q-e/-/releases/qe-6.5) (=6.5)
 - Compiler ([GNU](https://gcc.gnu.org/fortran/)(Recommended) or [Intel](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/fortran-compiler.html))

## Installation
### Pip
   Using pip can easy install the release version (serial) of QEpy from [PyPI](https://pypi.org/project/qepy).

```shell
python -m pip install qepy
```

### Source

 - **QE**

	All source codes should be compiled with the `-fPIC` compuiler option. Add `-fPIC` to the configuration options. E.g.,

     ```shell
	 ./configure CFLAGS=-fPIC FFLAGS=-fPIC try_foxflags=-fPIC MPIF90=mpif90
	  make pwall
	  export qedir=`pwd`
     ```

 - **QEpy**

     ```shell
	 git clone --recurse-submodules https://gitlab.com/shaoxc/qepy.git
     oldxml=yes ldau=yes tddft=yes python -m pip install -U ./qepy
	 ```

## Manual and Tutorials
  See [QEpy's website](http://qepy.rutgers.edu) for details.
