#!/usr/bin/env python3
import unittest
import numpy as np
from qepy.calculator import QEpyCalculator
import pathlib


class Test(unittest.TestCase):
    def test_md(self):
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
        except Exception:
            comm = None

        import ase.io
        from ase.io.trajectory import Trajectory
        from ase import units
        from ase.md.verlet import VelocityVerlet

        path = pathlib.Path(__file__).resolve().parent / 'DATA'
        inputfile = path / 'qe_in.in'

        calc = QEpyCalculator(comm = comm, inputfile = inputfile)
        atoms = ase.io.read(inputfile, format='espresso-in')
        atoms.calc = calc

        dyn = VelocityVerlet(atoms, 2 * units.fs)
        # traj = Trajectory("al_md_3.traj", "w", atoms)
        # dyn.attach(traj.write, interval=1)
        dyn.run(3)
        calc.stop()

        atoms_fin = ase.io.read(path / 'al_md_3.traj', index=-1)
        self.assertTrue(np.isclose(atoms.get_momenta(), atoms_fin.get_momenta(), atol=1.e-3).all())


if __name__ == "__main__":
    unittest.main()
