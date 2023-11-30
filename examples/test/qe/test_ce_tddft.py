import numpy as np
import qepy
import qepy_modules
import qepy_pw
import qepy_cetddft
import unittest
import pathlib
import shutil
import pytest

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
        qepy_pw.qepy_pwscf(inputfile, commf)
        qepy_pw.electrons()

        conv_flag = bool(qepy_modules.control_flags.get_conv_elec())
        self.assertTrue(conv_flag)

        etotal = qepy_pw.ener.get_etot()
        self.assertTrue(np.isclose(etotal, -552.93477389, atol = 1E-6))

        qepy_pw.punch('all')

    def test_1_tddft_continue(self):
        qepy_cetddft.qepy_tddft_readin(inputfile)
        qepy_cetddft.qepy_tddft_main_setup()
        qepy_cetddft.qepy_molecule_optical_absorption()
        qepy_pw.qepy_stop_run(0, print_flag=0, what='no', finalize=False)
        qepy_cetddft.qepy_stop_tddft(0)
        qepy.core.qepy_clean_saved()

    def test_2_tddft_iterative(self):
        embed = qepy_pw.qepy_common.embed_base()
        qepy_pw.qepy_common.set_embed(embed)
        qepy_cetddft.qepy_tddft_main_initial(inputfile, commf)
        embed.tddft.iterative = True
        qepy_pw.read_file()
        qepy_cetddft.qepy_tddft_main_setup()

        for i in range(5):
            qepy_cetddft.qepy_molecule_optical_absorption()
        dip = qepy_cetddft.qepy_tddft_common.get_array_dipole().copy()
        embed.tddft.finish = True
        qepy_cetddft.qepy_molecule_optical_absorption()
        qepy_cetddft.qepy_stop_tddft(0)

        assert abs(dip[0, 0] - 0.54355)<1E-3

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
