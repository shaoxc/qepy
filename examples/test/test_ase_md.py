from qepy.calculator import QEpyCalculator
import numpy as np
import pathlib

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

path = pathlib.Path(__file__).resolve().parent / 'DATA'
inputfile = path / 'qe_in.in'

def test_md():
    import ase.io
    from ase import units
    from ase.md.verlet import VelocityVerlet

    atoms = QEpyCalculator(inputfile = inputfile, comm = comm, from_file = True).atoms

    dyn = VelocityVerlet(atoms, 2 * units.fs)
    # from ase.io.trajectory import Trajectory
    # traj = Trajectory("al_md_3.traj", "w", atoms)
    # dyn.attach(traj.write, interval=1)
    dyn.run(3)
    atoms.calc.stop()

    atoms_fin = ase.io.read(path / 'al_md_3.traj', index=-1)
    assert np.allclose(atoms.get_momenta(), atoms_fin.get_momenta(), atol=1.e-3)
