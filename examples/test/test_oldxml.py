import qepy
from qepy.driver import Driver
import numpy as np
import pathlib
import pytest

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

@pytest.mark.skipif(not hasattr(qepy, 'oldxml_read_file'), reason="requires oldxml with qe-6.5")
def test_oldxml():
    if comm and comm.size > 1 : return
    path = pathlib.Path(__file__).resolve().parent / 'DATA/oldxml'
    driver = Driver(comm = comm, prefix = 'al_oldxml', outdir=path, task = 'nscf')
    energy = driver.get_energy()
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop(what = 'no')

@pytest.mark.skipif(not hasattr(qepy, 'oldxml_read_file'), reason="requires oldxml with qe-6.5")
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
