from qepy.driver import Driver
import numpy as np
import pathlib

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

def test_oldxml():
    path = pathlib.Path(__file__).resolve().parent / 'DATA/oldxml'
    driver = Driver(comm = comm, prefix = 'al_oldxml', outdir=path, task = 'nscf')
    energy = driver.get_energy()
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop(what = 'no')

def test_oldxml_collect():
    path = pathlib.Path(__file__).resolve().parent / 'DATA/oldxml_collect'
    driver = Driver(comm = comm, prefix = 'al_oldxml', outdir=path, task = 'nscf')
    wfc = path / 'al_oldxml.wfc1'
    if wfc.is_file(): wfc.unlink()
    wfc = path / 'al_oldxml.wfc'
    if wfc.is_file(): wfc.unlink()
    energy = driver.get_energy()
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop(what = 'no')


if __name__ == "__main__":
    test_oldxml()
    test_oldxml_collect()
