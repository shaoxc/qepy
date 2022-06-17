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
    @classmethod
    def setUpClass(cls):
        driver = QEpyDriver(inputfile, comm)
        driver.scf()
        converged = driver.check_convergence()
        cls.assertTrue(cls, converged)
        #
        energy = driver.get_energy()
        cls.assertTrue(cls, np.isclose(energy, -552.93477389, rtol = 1E-6))
        driver.save()

    def test_1_tddft_continue(self):
        driver = QEpyDriver(inputfile, comm, task = 'optical')
        driver.scf()
        dip = driver.get_dipole_tddft()
        driver.stop()
        assert(abs(dip[0, 0] - 0.56199)<1E-3)

    def test_2_tddft_iterative(self):
        driver = QEpyDriver(inputfile, comm, task = 'optical', iterative = True)
        for i in range(5):
            driver.diagonalize()
        dip = driver.get_dipole_tddft()
        driver.stop()
        assert(abs(dip[0, 0] - 0.54355)<1E-3)

    def test_3_tddft_restart(self):
        driver = QEpyDriver(inputfile, comm, task = 'optical', iterative = True)
        # restart from 5
        driver.tddft_restart(istep=5)
        for i in range(5):
            driver.diagonalize()
        dip = driver.get_dipole_tddft()
        driver.stop()
        assert(abs(dip[0, 0] - 0.56199)<1E-3)

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
