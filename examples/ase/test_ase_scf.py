#!/usr/bin/env python3
import qepy
import time

try:
    from mpi4py import MPI
except Exception:
    comm = None
else:
    comm = MPI.COMM_WORLD

from qepy.calculator import QEpyCalculator

inputfile = 'qe_in.in'

calc = QEpyCalculator(comm = comm, inputfile = inputfile)

get_potential_energy      = calc.get_potential_energy()
get_forces                = calc.get_forces()[0]
get_stress                = calc.get_stress()[0]
get_density               = calc.get_density()
get_bz_k_points           = calc.get_bz_k_points()[:, 0]
get_effective_potential   = calc.get_effective_potential()[0]
get_eigenvalues           = calc.get_eigenvalues()[0]
get_fermi_level           = calc.get_fermi_level()
get_ibz_k_points          = calc.get_ibz_k_points()[:, 0]
get_k_point_weights       = calc.get_k_point_weights()[0]
get_magnetic_moment       = calc.get_magnetic_moment()
get_number_of_bands       = calc.get_number_of_bands()
get_number_of_grid_points = calc.get_number_of_grid_points()
get_number_of_spins       = calc.get_number_of_spins()
get_occupation_numbers    = calc.get_occupation_numbers()[0]
get_pseudo_density        = calc.get_pseudo_density()[0]
get_pseudo_wave_function  = calc.get_pseudo_wave_function()[0, 0]
get_spin_polarized        = calc.get_spin_polarized()
get_xc_functional         = calc.get_xc_functional()

time.sleep(0.5)

if calc.rank == 0 :
    print('ncharge:', get_density.sum()*calc.atoms.get_volume()/get_density.shape[0], flush=True)
    print("calc.get_potential_energy()           =", get_potential_energy      ,  flush = True)
    print("calc.get_forces()[0]                  =", get_forces                ,  flush = True)
    print("calc.get_stress()[0]                  =", get_stress                ,  flush = True)
    print("calc.get_density()[0]                 =", get_density[0]            ,  flush = True)
    print("calc.get_bz_k_points()[:, 0]          =", get_bz_k_points           ,  flush = True)
    print("calc.get_effective_potential()[0]     =", get_effective_potential   ,  flush = True)
    print("calc.get_eigenvalues()[0]             =", get_eigenvalues           ,  flush = True)
    print("calc.get_fermi_level()                =", get_fermi_level           ,  flush = True)
    print("calc.get_ibz_k_points()[:, 0]         =", get_ibz_k_points          ,  flush = True)
    print("calc.get_k_point_weights()[0]         =", get_k_point_weights       ,  flush = True)
    print("calc.get_magnetic_moment()            =", get_magnetic_moment       ,  flush = True)
    print("calc.get_number_of_bands()            =", get_number_of_bands       ,  flush = True)
    print("calc.get_number_of_grid_points()      =", get_number_of_grid_points ,  flush = True)
    print("calc.get_number_of_spins()            =", get_number_of_spins       ,  flush = True)
    print("calc.get_occupation_numbers()[0]      =", get_occupation_numbers    ,  flush = True)
    print("calc.get_pseudo_density()[0]          =", get_pseudo_density        ,  flush = True)
    print("calc.get_pseudo_wave_function()[0, 0] =", get_pseudo_wave_function  ,  flush = True)
    print("calc.get_spin_polarized()             =", get_spin_polarized        ,  flush = True)
    print("calc.get_xc_functional()              =", get_xc_functional         ,  flush = True)
