from qepy.driver import Driver
import numpy as np
import pathlib

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

path = pathlib.Path(__file__).resolve().parent / 'DATA'
inputfile = path / 'qe_in.in'

def test_scf_iter():
    driver = Driver(inputfile, comm, iterative = True)
    for i in range(60):
        driver.diagonalize()
        driver.mix()
        converged = driver.check_convergence()
        if converged : break
    #
    energy = driver.get_energy()
    forces = driver.get_forces()
    stress = driver.get_stress()
    if driver.is_root :
        print('converged :\n', converged)
        print('energy :\n', energy)
        print('forces :\n', forces)
        print('stress :\n', stress)
    assert converged
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    assert np.isclose(forces[0, 0], -0.00835135, atol = 1E-3)
    assert np.isclose(stress[1, 1], -0.00256059, atol = 1E-3)
    driver.stop()
