import qepy
import qepy_pw
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
inputfile = path / 'qe_fake.in'


class Test(unittest.TestCase):
    def test_pw_comm(self):
        # return # Take some time to finish
        for i in range(1000): # Error happened ~450
            qepy_pw.qepy_pwscf(inputfile, commf)
            # qepy_pw.qepy_mod.qepy_set_stdout('/dev/null')
            qepy_pw.qepy_mod.qepy_set_stdout('al.out')
            qepy_pw.qepy_stop_run(0, what = 'no')
            qepy_pw.qepy_mod.qepy_write_stdout(f'Test step = {i}')
            print('istep = ', i, flush = True)

    def test_initial_comm(self):
        return # read_file take some time
        inputobj = qepy_pw.qepy_common.input_base()
        inputobj.prefix = 'tmp'
        if commf : inputobj.my_world_comm = commf
        for i in range(50000): # Error happened ~16000
            # qepy_pw.qepy_mod.qepy_set_stdout('/dev/null')
            qepy_pw.qepy_mod.qepy_set_stdout('al.out')
            qepy_pw.qepy_initial(inputobj)
            # qepy_pw.read_file()
            qepy_pw.qepy_stop_run(0, what = 'no')
            qepy_pw.qepy_mod.qepy_write_stdout(f'Test step = {i}')

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
