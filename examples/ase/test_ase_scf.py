#!/usr/bin/env python3
import qepy

try:
    from mpi4py import MPI
except Exception:
    comm = None
else:
    comm = MPI.COMM_WORLD

from qepy.calculator import QEpyCalculator

inputfile = 'qe_in.in'

calc = QEpyCalculator(comm = comm, inputfile = inputfile)

energy = calc.get_potential_energy()
forces = calc.get_forces()
stress = calc.get_stress()
density = calc.get_density()

if calc.rank == 0 :
    print('qepy -> energy :', energy)
    print('qepy -> forces :\n', forces)
    print('qepy -> stress :\n', stress)
    print('qepy -> density :', density.sum()*calc.atoms.get_volume()/density.size)
