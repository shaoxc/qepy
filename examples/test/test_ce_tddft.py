from qepy.driver import Driver
import numpy as np
import pathlib
import pytest

path = pathlib.Path(__file__).resolve().parent / 'DATA'
inputfile = path / 'qe_in.in'

def test_0_scf():
    pytest.importorskip("qepy_cetddft")
    # Run scf
    driver = Driver(inputfile, comm=True)
    driver.scf()
    converged = driver.check_convergence()
    energy = driver.get_energy()
    driver.save()
    #
    energy = driver.get_energy()
    if driver.is_root :
        print('converged :\n', converged)
        print('energy :\n', energy)
    assert converged
    assert np.isclose(energy, -552.93477389, atol = 1E-6)

def test_1_tddft_continue():
    pytest.importorskip("qepy_cetddft")
    # Run TDDFT after scf, without stop
    driver = Driver(inputfile, comm=True, task = 'optical', progress = True)
    driver.scf()
    dipole = driver.get_dipole_tddft()
    if driver.is_root :
        print('dipople:\n', dipole)
    assert np.isclose(dipole[0, 0], 0.56199, atol = 1E-3)
    driver.stop()

def test_2_tddft_iterative():
    pytest.importorskip("qepy_cetddft")
    # Run TDDFT
    driver = Driver(inputfile, comm=True, task = 'optical', iterative = True)
    for i in range(5):
        driver.diagonalize()
    dipole = driver.get_dipole_tddft()
    if driver.is_root :
        print('dipople:\n', dipole)
    assert np.isclose(dipole[0, 0], 0.54355, atol = 1E-3)
    driver.stop()

def test_3_tddft_restart():
    pytest.importorskip("qepy_cetddft")
    driver = Driver(inputfile, comm=True, task = 'optical', iterative = True)
    # restart from 5
    driver.tddft_restart(istep=5)
    for i in range(5):
        driver.diagonalize()
    dipole = driver.get_dipole_tddft()
    if driver.is_root :
        print('dipople:\n', dipole)
    assert np.isclose(dipole[0, 0], 0.56199, atol = 1E-3)
    driver.stop()


if __name__ == "__main__":
    tests = [item for item in globals() if item.startswith('test_')]
    for func in sorted(tests):
        globals()[func]()
