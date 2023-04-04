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

        qepy.qepy_stop_run(0, what = 'all')

    def test_1_read(self):
        inputobj = qepy.qepy_common.input_base()
        inputobj.prefix = 'tmp'
        if commf : inputobj.my_world_comm = commf

        qepy.qepy_initial(inputobj)

        qepy.read_file()

        embed = qepy.qepy_common.embed_base()
        qepy.qepy_common.set_embed(embed)
        qepy.qepy_calc_energies()
        self.assertTrue(np.isclose(embed.etotal, -552.93477389, atol = 1E-6))
        qepy.qepy_stop_run(0, what = 'no')

    def test_2_read_pw(self):
        qepy.qepy_pwscf(inputfile, commf)
        embed = qepy.qepy_common.embed_base()
        qepy.qepy_common.set_embed(embed)

        qepy.qepy_mod.qepy_restart_from_xml()
        if qepy.basis.get_starting_pot().strip() != 'file' :
            qepy.basis.set_starting_pot('file')
            qepy.potinit()
        if qepy.basis.get_starting_wfc().strip() != 'file' :
            qepy.basis.set_starting_wfc('file')
            qepy.wfcinit()

        qepy.qepy_calc_energies()
        self.assertTrue(np.isclose(embed.etotal, -552.93477389, atol = 1E-6))
        qepy.qepy_stop_run(0, what = 'no')

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
