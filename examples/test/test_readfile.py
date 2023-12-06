from qepy.driver import Driver
import numpy as np
import pathlib

path = pathlib.Path(__file__).resolve().parent / 'DATA'
inputfile = path / 'qe_in.in'

def test_0_scf():
    driver = Driver(inputfile, comm=True)
    driver.scf()
    converged = driver.check_convergence()
    energy = driver.get_energy()
    assert converged
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop()

def test_1_read():
    driver = Driver(prefix = 'tmp', task = 'nscf', comm=True)
    energy = driver.get_energy()
    if driver.is_root :
        print('energy :\n', energy)
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop(what = 'no')

def test_2_read_pw():
    driver = Driver(inputfile, comm=True)
    driver.pwscf_restart()
    energy = driver.get_energy()
    if driver.is_root :
        print('energy :\n', energy)
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    driver.stop(what = 'no')


if __name__ == "__main__":
    tests = [item for item in globals() if item.startswith('test_')]
    for func in sorted(tests):
        globals()[func]()
