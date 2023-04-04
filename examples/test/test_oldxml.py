import qepy
from qepy.driver import Driver
import numpy as np
import pathlib

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

def test_oldxml():
    if comm and comm.size > 1 : return
    path = pathlib.Path(__file__).resolve().parent / 'DATA/oldxml'
    driver = Driver(comm = comm, prefix = 'al_oldxml', outdir=path, task = 'nscf')
    energy = driver.get_energy()
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop(what = 'no')

def test_oldxml_collect():
    path = pathlib.Path(__file__).resolve().parent / 'DATA/oldxml_collect'
    driver = Driver(comm = comm, prefix = 'al_oldxml', outdir=path, task = 'nscf')
    energy = driver.get_energy()
    #
    if driver.is_root :
        for f in path.glob('al_oldxml.wfc*'):
            f.unlink()
    #
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop(what = 'no')


if __name__ == "__main__":
    test_oldxml()
    test_oldxml_collect()
