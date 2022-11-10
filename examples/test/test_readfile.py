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
    driver = Driver(inputfile, comm)
    driver.scf()
    converged = driver.check_convergence()
    energy = driver.get_energy()
    assert converged
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop()

def test_1_read():
    driver = Driver(comm = comm, prefix = 'tmp', task = 'nscf')
    energy = driver.get_energy()
    if driver.is_root :
        print('energy :\n', energy)
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop(what = 'no')

def test_2_read_pw():
    driver = Driver(inputfile, comm)
    driver.pwscf_restart()
    energy = driver.get_energy()
    if driver.is_root :
        print('energy :\n', energy)
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop(what = 'no')
