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
    def test_0_scf(self):
        qepy.qepy_pwscf(inputfile, commf)
        qepy.electrons()

        conv_flag = bool(qepy.control_flags.get_conv_elec())
        self.assertTrue(conv_flag)

        etotal = qepy.ener.get_etot()
        self.assertTrue(np.isclose(etotal, -552.93477389, atol = 1E-6))

        qepy.punch('all')

    def test_1_tddft_continue(self):
        qepy.qepy_tddft_readin(inputfile)
        qepy.qepy_tddft_main_setup()
        qepy.qepy_molecule_optical_absorption()
        qepy.qepy_stop_run(0, print_flag=0, what='no', finalize=False)
        qepy.qepy_stop_tddft(0)

    def test_2_tddft_iterative(self):
        embed = qepy.qepy_common.embed_base()
        qepy.qepy_tddft_main_initial(inputfile, commf, embed = embed)
        qepy.read_file()
        embed.tddft.iterative = True
        qepy.qepy_tddft_main_setup()

        for i in range(5):
            qepy.qepy_molecule_optical_absorption()
        dip = qepy.qepy_tddft_common.get_array_dipole().copy()
        embed.tddft.finish = True
        qepy.qepy_molecule_optical_absorption()
        qepy.qepy_stop_tddft(0)

        assert(abs(dip[0, 0] - 0.54355)<1E-3)

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
