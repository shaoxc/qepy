#!/usr/bin/env python3
import qepy
from qepy.calculator import QEpyCalculator, units

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

import ase.io
from ase.io.trajectory import Trajectory

from ase.optimize import BFGS, LBFGS, FIRE
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG
from ase.constraints import StrainFilter, UnitCellFilter
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry

inputfile = 'vcrelax.in'

calc = QEpyCalculator(comm = comm, inputfile = inputfile, lmovecell = True)

atoms = ase.io.read(inputfile, format='espresso-in')

atoms.set_calculator(calc)
atoms.set_constraint(FixSymmetry(atoms))

# af = atoms 

pressure = 10 # GPa
af = UnitCellFilter(atoms, scalar_pressure = pressure*units['J']/10**21)

trajfile = 'opt.traj'
opt = LBFGS(af, trajectory = trajfile, memory = 8, use_line_search = False)

opt.run(fmax = 0.01)

if comm and comm.rank == 0 :
    traj = Trajectory(trajfile)
    ase.io.write('opt.vasp', traj[-1], direct = True, long_format=True, vasp5 = True)
    ase.io.write('opt.in', traj[-1], format = 'espresso-in')
