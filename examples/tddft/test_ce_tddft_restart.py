from qepy.driver import Driver
import numpy as np
import pathlib

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

inputfile = 'qe_in.in'

def main():
    # Run TDDFT iterative
    driver = Driver(inputfile, comm, task = 'optical', iterative = True)
    for i in range(5):
        driver.diagonalize()
    dipole = driver.get_dipole_tddft()
    if driver.is_root :
        print('dipople:\n', dipole)
    driver.stop()
    # Restart the TDDFT
    driver = Driver(inputfile, comm, task = 'optical', iterative = True)
    # restart from 5
    driver.tddft_restart(istep=5)
    for i in range(5):
        driver.diagonalize()
    dipole = driver.get_dipole_tddft()
    if driver.is_root :
        print('dipople:\n', dipole)
    driver.stop()


if __name__ == "__main__":
    main()
