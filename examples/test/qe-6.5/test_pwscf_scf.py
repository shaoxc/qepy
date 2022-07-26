import numpy as np
import qepy
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
        qepy.qepy_pwscf(inputfile, commf)
        embed = qepy.qepy_common.embed_base()
        # embed.ldescf = True # add scf correction energy
        embed.iterative = True
        for i in range(60):
            embed.mix_coef = -1.0
            qepy.qepy_electrons_scf(2, 0, embed)
            embed.mix_coef = 0.7
            qepy.qepy_electrons_scf(2, 0, embed)
            if qepy.control_flags.get_conv_elec() : break

        conv_flag = bool(qepy.control_flags.get_conv_elec())
        self.assertTrue(conv_flag)

        etotal = embed.etotal
        self.assertTrue(np.isclose(etotal, -552.93477389, atol = 1E-6))

        qepy.qepy_stop_run(0, what = 'no')

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
