import numpy as np
import qepy
from qepy import qepy_pw
from qepy import qepy_modules
import unittest
import pathlib
import shutil
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    commf = comm.py2f()
except Exception:
    comm = None
    commf = None

path = pathlib.Path(__file__).resolve().parent / 'DATA'
inputfile = path / 'qe_in.in'


class Test(unittest.TestCase):
    def test_scf(self):
        embed = qepy_pw.qepy_common.embed_base()
        qepy_pw.qepy_pwscf(inputfile, commf, embed = embed)
        qepy_pw.qepy_electrons_scf(2, 0)

        conv_flag = bool(qepy_modules.control_flags.get_conv_elec())
        self.assertTrue(conv_flag)

        etotal = embed.etotal
        self.assertTrue(np.isclose(etotal, -552.93477389, atol = 1E-6))

        qepy_pw.qepy_forces(0)
        forces = qepy_pw.force_mod.get_array_force().T
        self.assertTrue(np.isclose(forces[0, 0], -0.00835135, atol = 1E-3))

        stress = np.ones((3, 3), order='F')
        qepy_pw.qepy_stress(stress)
        self.assertTrue(np.isclose(stress[1, 1], -0.00256059, atol = 1E-3))

        qepy_pw.qepy_stop_run(0, what = 'no')

    def tearDown(self):
        if comm and comm.rank == 0 :
            path = pathlib.Path('.')
            for f in path.glob('al.*'):
                if f.is_file():
                    f.unlink()
                else :
                    shutil.rmtree(f)


if __name__ == "__main__":
    unittest.main()
