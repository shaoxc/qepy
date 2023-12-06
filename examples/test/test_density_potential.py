#!/usr/bin/env python
from qepy.driver import Driver
import numpy as np

qe_options = {
    '&control': {
        'calculation': "'scf'",
        'pseudo_dir': "'./DATA/'",
        'prefix' : "'tmp'"
    },
    '&system': {
        'ibrav' : 0,
        'degauss': 0.005,
        'ecutwfc': 30,
        'nat': 1,
        'ntyp': 1,
        'occupations': "'smearing'"
    },
    '&electrons': {
        'conv_thr' : 1e-10
    },
    'atomic_positions crystal': ['Al    0.0  0.0  0.0'],
    'atomic_species': ['Al  26.98 Al_ONCV_PBE-1.2.upf'],
    'k_points automatic': ['2 2 2 0 0 0'],
    'cell_parameters angstrom':[
        '0.     2.025  2.025',
        '2.025  0.     2.025',
        '2.025  2.025  0.   '],
}

etotal = -137.929763
saved = {}

def test_0_scf_potential():
    driver=Driver(qe_options=qe_options, comm=True, logfile='tmp.0.out')
    driver.scf()
    rho = driver.get_density().copy()
    energy = driver.get_energy()
    xc = driver.get_exchange_correlation()
    if driver.is_root :
        print('total:', energy)
        print('exc:', xc[1])
    assert np.isclose(energy, etotal, atol = 1E-5)
    assert np.isclose(xc[1], -15.5051464238815, atol=1E-5)
    driver.set_external_potential(xc[0], exttype=('xc'))
    driver.scf()
    extotal = driver.calc_energy()
    assert np.isclose(extotal, -122.424617536980, atol = 1E-5)
    saved['rho'] = rho

def test_1_density_nscf():
    driver=Driver(qe_options=qe_options, comm=True, logfile='tmp.1.out')
    #
    driver.set_density(saved['rho'])
    driver.qepy_pw.qepy_mod.qepy_calc_effective_potential()
    driver.qepy_modules.control_flags.set_ethr(1E-9)
    driver.non_scf()
    #
    energy = driver.calc_energy()
    if driver.is_root :
        print('etotal:', energy)
    assert np.isclose(energy, etotal, atol = 1E-5)


if __name__ == "__main__":
    tests = [item for item in globals() if item.startswith('test_')]
    for func in sorted(tests):
        globals()[func]()
