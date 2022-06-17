from qepy.driver import Driver
import numpy as np
import pathlib

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

path = pathlib.Path(__file__).resolve().parent / 'DATA'
inputfile = path / 'qe_in.in'

def test_0_scf():
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
    #
    assert converged
    assert np.isclose(energy, -552.93477389, rtol = 1E-6)

def test_1_tddft_continue():
    # Run TDDFT
    driver = Driver(inputfile, comm, task = 'optical')
    driver.scf()
    dipole = driver.get_dipole_tddft()
    if driver.is_root :
        print('dipople:\n', dipole)
    driver.stop()
    #
    assert(abs(dipole[0, 0] - 0.56199)<1E-3)

def test_2_tddft_iterative():
    driver = Driver(inputfile, comm, task = 'optical', iterative = True)
    for i in range(5):
        driver.diagonalize()
    dipole = driver.get_dipole_tddft()
    if driver.is_root :
        print('dipople:\n', dipole)
    driver.stop()
    #
    assert(abs(dipole[0, 0] - 0.54355)<1E-3)

def test_3_tddft_restart():
    driver = Driver(inputfile, comm, task = 'optical', iterative = True)
    # restart from 5
    driver.tddft_restart(istep=5)
    for i in range(5):
        driver.diagonalize()
    dipole = driver.get_dipole_tddft()
    if driver.is_root :
        print('dipople:\n', dipole)
    driver.stop()
    #
    assert(abs(dipole[0, 0] - 0.56199)<1E-3)
