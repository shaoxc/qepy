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
    # Run scf
    driver = Driver(inputfile, comm)
    driver.scf()
    converged = driver.check_convergence()
    energy = driver.get_energy()
    driver.save()
    #
    energy = driver.get_energy()
    if driver.is_root :
        print('converged :\n', converged)
        print('energy :\n', energy)
    # Run TDDFT
    driver = Driver(inputfile, comm, task = 'optical')
    driver.scf()
    dipole = driver.get_dipole_tddft()
    if driver.is_root :
        print('dipople:\n', dipole)
    driver.stop()


if __name__ == "__main__":
    main()
