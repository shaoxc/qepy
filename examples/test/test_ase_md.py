#!/usr/bin/env python3
import unittest
import numpy as np
from qepy.calculator import QEpyCalculator
import pathlib
import shutil
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

class Test(unittest.TestCase):
    def test_md(self):
        import ase.io
        from ase import units
        from ase.md.verlet import VelocityVerlet

        path = pathlib.Path(__file__).resolve().parent / 'DATA'
        inputfile = path / 'qe_in.in'

        calc = QEpyCalculator(comm = comm, inputfile = inputfile)
        atoms = ase.io.read(inputfile, format='espresso-in')
        atoms.calc = calc

        dyn = VelocityVerlet(atoms, 2 * units.fs)
        # from ase.io.trajectory import Trajectory
        # traj = Trajectory("al_md_3.traj", "w", atoms)
        # dyn.attach(traj.write, interval=1)
        dyn.run(3)
        calc.stop()

        atoms_fin = ase.io.read(path / 'al_md_3.traj', index=-1)
        self.assertTrue(np.isclose(atoms.get_momenta(), atoms_fin.get_momenta(), atol=1.e-3).all())

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
