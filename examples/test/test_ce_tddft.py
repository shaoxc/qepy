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

class Test(unittest.TestCase):
    def test_0_scf(self):
        global commf
        path = pathlib.Path(__file__).resolve().parent / 'DATA'
        fname = path / 'qe_in.in'

        qepy.qepy_pwscf(fname, commf)
        qepy.electrons()

        conv_flag = bool(qepy.control_flags.get_conv_elec())
        self.assertTrue(conv_flag)

        etotal = qepy.ener.get_etot()
        self.assertTrue(np.isclose(etotal, -552.93477389, rtol = 1E-6))

        qepy.punch('all')

    def test_1_tddft_continue(self):
        global commf
        path = pathlib.Path(__file__).resolve().parent / 'DATA'
        fname = path / 'qe_in.in'

        embed = qepy.qepy_common.embed_base()
        qepy.qepy_tddft_readin(fname)
        qepy.qepy_tddft_main_setup(embed)
        qepy.qepy_molecule_optical_absorption(embed)
        qepy.qepy_stop_run(0, what = 'no')
        qepy.qepy_stop_tddft(0)

    def test_2_tddft_iterative(self):
        global commf
        path = pathlib.Path(__file__).resolve().parent / 'DATA'
        fname = path / 'qe_in.in'

        qepy.qepy_tddft_main_initial(fname, commf)
        qepy.read_file()
        embed = qepy.qepy_common.embed_base()
        embed.tddft.iterative = True
        qepy.qepy_tddft_main_setup(embed)

        for i in range(5):
            qepy.qepy_molecule_optical_absorption(embed)
        dip = qepy.qepy_tddft_common.get_array_dipole().copy()
        embed.tddft.finish = True
        qepy.qepy_molecule_optical_absorption(embed)
        qepy.qepy_stop_tddft(0)

        assert(abs(dip[0, 0] - 0.5435)<1E-3)

    def test_9_clean(self):
        if comm and comm.rank == 0 :
            path = pathlib.Path('.')
            for f in path.glob('al.*'):
                if f.is_file():
                    f.unlink()
                else :
                    shutil.rmtree(f)


if __name__ == "__main__":
    unittest.main()
