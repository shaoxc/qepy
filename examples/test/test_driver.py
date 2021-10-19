import numpy as np
import qepy
from qepy.driver import QEpyDriver
import unittest
import pathlib
import shutil

class Test(unittest.TestCase):
    def test_driver(self):
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            comm = comm.py2f()
        except Exception:
            comm = None

        path = pathlib.Path(__file__).resolve().parent / 'DATA'
        fname = path / 'qe_in.in'

        driver = QEpyDriver(fname, comm)

        for i in range(60):
            driver.diagonalize()
            driver.mix(mix_coef = 0.7)
            if driver.check_convergence(): break

        energy = driver.get_energy()
        self.assertTrue(np.isclose(energy, -552.93477389, rtol = 1E-6))

        forces = driver.get_forces()
        self.assertTrue(np.isclose(forces[0, 0], -0.00835135, rtol = 1E-3))

        stress = driver.get_stress()
        self.assertTrue(np.isclose(stress[1, 1], -0.00256059, rtol = 1E-3))

        driver.stop()

    def tearDown(self):
        path = pathlib.Path('.')
        for f in path.glob('al.*'):
            if f.is_file():
                f.unlink()
            else :
                shutil.rmtree(f)


if __name__ == "__main__":
    unittest.main()
