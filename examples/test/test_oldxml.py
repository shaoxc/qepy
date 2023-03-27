from qepy.driver import Driver
import numpy as np
import pathlib

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

path = pathlib.Path(__file__).resolve().parent / 'DATA/oldxml'

def test_oldxml():
    driver = Driver(comm = comm, prefix = 'tmp', outdir=path, task = 'nscf')
    energy = driver.get_energy()
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop()
