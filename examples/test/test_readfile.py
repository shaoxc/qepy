import numpy as np
import qepy
import unittest
import pathlib
import shutil

class Test(unittest.TestCase):
    def test_0_scf(self):
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            comm = comm.py2f()
        except Exception:
            comm = None
        path = pathlib.Path(__file__).resolve().parent / 'DATA'
        fname = path / 'qe_in.in'
        qepy.qepy_pwscf(fname, comm)
        qepy.electrons()

        conv_flag = bool(qepy.control_flags.get_conv_elec())
        self.assertTrue(conv_flag)

        etotal = qepy.ener.get_etot()
        self.assertTrue(np.isclose(etotal, -552.93477389, rtol = 1E-6))

        qepy.qepy_stop_run(0, what = 'all')

    def test_1_read(self):
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            comm = comm.py2f()
        except Exception:
            comm = None

        inputobj = qepy.qepy_common.input_base()
        inputobj.prefix = 'al'
        if comm : inputobj.my_world_comm = comm

        qepy.qepy_initial(inputobj)

        qepy.read_file()

        embed = qepy.qepy_common.embed_base()
        qepy.qepy_calc_energies(embed)
        self.assertTrue(np.isclose(embed.etotal, -552.93477389, rtol = 1E-6))
        qepy.qepy_stop_run(0, what = 'no')

    def test_2_read_pw(self):
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            comm = comm.py2f()
        except Exception:
            comm = None

        path = pathlib.Path(__file__).resolve().parent / 'DATA'
        fname = path / 'qe_in.in'

        qepy.qepy_pwscf(fname, comm)
        embed = qepy.qepy_common.embed_base()

        qepy.qepy_pw_restart_new.qepy_read_xml_file(alloc=False)
        if qepy.basis.get_starting_pot().strip() != 'file' :
            qepy.qepy_potinit(starting = 'file')
        if qepy.basis.get_starting_wfc().strip() != 'file' :
            qepy.qepy_wfcinit(starting = 'file')

        qepy.qepy_calc_energies(embed)
        self.assertTrue(np.isclose(embed.etotal, -552.93477389, rtol = 1E-6))
        qepy.qepy_stop_run(0, what = 'no')

    def test_9_clean(self):
        path = pathlib.Path('.')
        for f in path.glob('al.*'):
            if f.is_file():
                f.unlink()
            else :
                shutil.rmtree(f)


if __name__ == "__main__":
    unittest.main()
