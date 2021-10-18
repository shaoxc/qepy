import numpy as np
import qepy
import unittest
import pathlib
import shutil

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    comm = comm.py2f()
except Exception:
    comm = None

class Test(unittest.TestCase):
    def test_0_scf(self):
        global comm
        path = pathlib.Path(__file__).resolve().parent / 'DATA'
        fname = path / 'qe_in.in'

        qepy.qepy_pwscf(fname, comm)
        qepy.electrons()

        conv_flag = bool(qepy.control_flags.get_conv_elec())
        self.assertTrue(conv_flag)

        etotal = qepy.ener.get_etot()
        self.assertTrue(np.isclose(etotal, -552.93477389, rtol = 1E-6))

        qepy.punch('all')

    def test_1_tddft_continue(self):
        global comm
        path = pathlib.Path(__file__).resolve().parent / 'DATA'
        fname = path / 'qe_in.in'

        embed = qepy.qepy_common.embed_base()
        qepy.qepy_tddft_readin(fname)
        qepy.qepy_tddft_main_setup(embed)
        qepy.qepy_molecule_optical_absorption(embed)
        qepy.qepy_stop_run(0, what = 'no')
        qepy.qepy_stop_tddft(0)

    def test_2_tddft_iterative(self):
        global comm
        path = pathlib.Path(__file__).resolve().parent / 'DATA'
        fname = path / 'qe_in.in'

        qepy.qepy_tddft_main_initial(fname, comm)
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

        assert(abs(dip[0, 0] - 0.54355185)<1E-6)

    def test_9_clean(self):
        path = pathlib.Path('.')
        for f in path.glob('al.*'):
            if f.is_file():
                f.unlink()
            else :
                shutil.rmtree(f)


if __name__ == "__main__":
    unittest.main()
