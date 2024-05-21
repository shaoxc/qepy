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

tol = 1e-8
etotal = -137.929763

def test_1_linear():
    driver=Driver(qe_options=qe_options, comm=True, iterative = True, logfile='tmp.1.out')
    rho = driver.get_density().copy()
    coef = 0.7
    for i in range(60):
        driver.diagonalize()
        rho_new = driver.get_density().copy()
        drho = rho_new-rho
        nc = np.abs(drho).sum()*driver.get_volume()/drho.size
        if driver.is_root :
            print('Iter: ',i, 'dN:', nc)
        if driver.nproc > 1: nc = driver.comm.bcast(nc)
        if nc < tol :
            driver.end_scf()
            break
        rho = (1-coef)*rho + coef * rho_new
        driver.set_density(rho)
    energy = driver.get_energy()
    if driver.is_root :
        print(energy)
    assert np.isclose(energy, etotal, atol = 1E-5)

def test_2_mix_potential():
    driver=Driver(qe_options=qe_options, comm=True, iterative = True, logfile='tmp.2.out')
    rho_new = driver.get_density().copy()
    coef = 0.7
    hartree = driver.get_hartree()
    xc = driver.get_exchange_correlation()
    v_hxc = hartree[0] + xc[0]
    for i in range(60):
        driver.diagonalize()
        rho_new, rho = driver.get_density().copy(), rho_new
        drho = rho_new-rho
        nc = np.abs(drho).sum()*driver.get_volume()/drho.size
        if driver.is_root :
            print('Iter: ',i, 'dN:', nc)
        if driver.nproc > 1: nc = driver.comm.bcast(nc)
        if nc < tol:
            driver.end_scf()
            break
        hartree = driver.get_hartree()
        xc = driver.get_exchange_correlation()
        v_hxc_new = hartree[0] + xc[0]
        v_hxc = (1-coef)*v_hxc + coef*v_hxc_new
        extene = hartree[1] + xc[1]
        driver.set_external_potential(v_hxc, exttype=('hartree', 'xc'), extene=extene)
    energy = driver.get_energy()
    if driver.is_root :
        print(energy)
    assert np.isclose(energy, etotal, atol = 1E-5)


if __name__ == "__main__":
    tests = [item for item in globals() if item.startswith('test_')]
    for func in sorted(tests):
        globals()[func]()
