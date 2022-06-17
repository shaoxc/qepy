import numpy as np
from qepy.driver import QEpyDriver
import unittest
import pathlib
import shutil

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

path = pathlib.Path(__file__).resolve().parent / 'DATA'
inputfile = path / 'qe_in.in'

class Test(unittest.TestCase):
    def test_0_scf(self):
        driver = QEpyDriver(inputfile, comm)
        driver.scf()
        converged = driver.check_convergence()
        self.assertTrue(converged)
        #
        energy = driver.get_energy()
        self.assertTrue(np.isclose(energy, -552.93477389, rtol = 1E-6))
        driver.stop()

    def test_1_read(self):
        driver = QEpyDriver(comm = comm, prefix = 'al', task = 'nscf')
        #
        energy = driver.get_energy()
        self.assertTrue(np.isclose(energy, -552.93477389, rtol = 1E-6))
        driver.stop(what = 'no')

    def test_2_read_pw(self):
        driver = QEpyDriver(inputfile, comm)
        driver.pwscf_restart()
        #
        energy = driver.get_energy()
        self.assertTrue(np.isclose(energy, -552.93477389, rtol = 1E-6))
        driver.stop(what = 'no')

    @classmethod
    def tearDownClass(cls):
        if comm and comm.rank == 0 :
            path = pathlib.Path('.')
            for f in path.glob('al.*'):
                if f.is_file():
                    f.unlink()
                else :
                    shutil.rmtree(f)


if __name__ == "__main__":
    unittest.main()
