import numpy as np
import qepy
from qepy import qepy_modules
from qepy import qepy_pw
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
        qepy_pw.qepy_common.set_embed(embed)
        qepy_pw.qepy_pwscf(inputfile, commf)
        # embed.ldescf = True # add scf correction energy
        embed.iterative = True
        for i in range(60):
            embed.mix_coef = -1.0
            qepy_pw.qepy_electrons_scf(2, 0)
            embed.mix_coef = 0.7
            qepy_pw.qepy_electrons_scf(2, 0)
            if qepy_modules.control_flags.get_conv_elec() : break

        conv_flag = bool(qepy_modules.control_flags.get_conv_elec())
        self.assertTrue(conv_flag)

        etotal = embed.etotal
        self.assertTrue(np.isclose(etotal, -552.93477389, atol = 1E-6))

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
