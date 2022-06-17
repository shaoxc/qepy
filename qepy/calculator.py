import numpy as np
from qepy.driver import QEpyDriver
from qepy import constants
from ase.calculators.calculator import Calculator
from ase.units import create_units
from ase.geometry import wrap_positions
import ase.io

units = create_units('2006').copy()
units['Bohr'] = constants.BOHR_RADIUS_SI * 1E10
units['Ry'] = constants.RYTOEV

class QEpyCalculator(Calculator):
    """QEpy calculator for ase"""
    implemented_properties=['energy', 'forces', 'stress']

    def __init__(self, atoms = None, inputfile = None, input_data = None,
            from_file = False, wrap = False, lmovecell = False,
            comm = None, task = 'scf', embed = None, ldescf = None, iterative = False, **kwargs):
        Calculator.__init__(self, atoms = atoms, **kwargs)
        self.optimizer = None
        self.atoms = atoms
        self.inputfile = inputfile
        self.wrap = wrap
        self.lmovecell = lmovecell
        self.from_file = from_file
        self.input_data = input_data
        #
        if self.inputfile is not None and self.atoms is None :
            self.atoms = ase.io.read(self.inputfile, format='espresso-in')
        if self.atoms is not None :
            self.atoms.calc = self
        #
        if self.from_file :
            self.basefile = None
        else :
            self.basefile = self.inputfile
            self.inputfile = 'qepy_input_tmp.in'
        #
        self.atoms_save = None
        self.driver = None
        self.iter = 0
        #
        self.qepy_options = {
                'comm' : comm,
                'ldescf' : ldescf,
                'iterative' : iterative,
                'task' : task,
                'embed' : embed,
                }
        self.qepy_options.update(kwargs)
        self.restart()

    def restart(self):
        self._energy = None
        self._forces = None
        self._stress = None

    def driver_initialise(self, atoms = None, **kwargs):
        atoms = atoms or self.atoms
        if self.basefile :
            if self.wrap : self.atoms.wrap()
            write2qe(self.inputfile, atoms, basefile = self.basefile, **self.parameters)
        self.driver = QEpyDriver(self.inputfile, **self.qepy_options)

    def update_atoms(self, atoms = None, **kwargs):
        atoms = atoms or self.atoms
        if self.wrap :
            positions = wrap_positions(atoms.positions, atoms.cell)
        else :
            positions = atoms.positions

        positions = positions /units['Bohr']
        if self.lmovecell :
            lattice = atoms.cell /units['Bohr']
        else :
            lattice = None
        self.driver.update_ions(positions=positions, lattice=lattice, **kwargs)

    def check_restart(self, atoms = None):
        restart = True

        if atoms is None or (self.atoms_save and atoms == self.atoms_save):
            restart = False

        if self.iter == 0 :
            restart = True

        if restart :
            if self.atoms is None :
                self.atoms = atoms
                self.atoms.calc = self
            if atoms is not None :
                self.atoms_save = atoms.copy()
            else :
                atoms = self.atoms
            self.restart()
            if self.iter > 0 :
                self.update_atoms(atoms)
            else :
                self.driver_initialise(atoms = atoms)

        return restart

    def update_optimizer(self, atoms = None, **kwargs):
        if self.check_restart(atoms):
            restart = True
            self.iter += 1
            self.driver.scf(**kwargs)
        else :
            restart = False
        return restart

    def stop(self, what = 'all', **kwargs):
        return self.driver.stop(what = 'all', **kwargs)

    @property
    def is_root(self):
        return self.driver.is_root

    @staticmethod
    def get_atoms_from_qepy():
        from ase.atoms import Atoms
        symbols = QEpyDriver.get_ions_symbols()
        positions = QEpyDriver.get_ions_positions() * units['Bohr']
        lattice = QEpyDriver.get_ions_lattice() * units['Bohr']
        atoms = Atoms(symbols = symbols, positions = positions, cell = lattice)
        return atoms

    def get_density(self, gather = True):
        """Return density array in real space."""
        return self.driver.get_density(gather=gather) / (units['Bohr'] ** 3)

    def get_wave_function(self, band=None, kpt=0):
        """Return wave-function array in real space."""
        return self.driver.get_wave_function(band=band, kpt=kpt)

    def get_number_of_k_points(self):
        return self.driver.get_number_of_k_points()

    def get_potential_energy(self, atoms = None, **kwargs):
        self.update_optimizer(atoms)
        return self.driver.get_energy(**kwargs) * units['Ry']

    def get_forces(self, atoms = None, icalc = 0):
        if self.update_optimizer(atoms) or self._forces is None:
            self._forces = self.driver.get_forces(icalc=icalc)* units['Ry'] / units['Bohr']
        return self._forces

    def get_stress(self, atoms = None):
        if self.update_optimizer(atoms) or self._stress is None:
            stress = self.driver.get_stress()
            self._stress = np.array(
                [stress[0, 0], stress[1, 1], stress[2, 2],
                 stress[1, 2], stress[0, 2], stress[0, 1]]) * -1 * units['Ry'] / (units['Bohr'] ** 3)
        return self._stress

    def get_number_of_bands(self):
        """Return the number of bands."""
        return self.driver.get_number_of_bands()

    def get_xc_functional(self):
        """Return the XC-functional identifier.

        'LDA', 'PBE', ..."""
        return self.driver.get_xc_functional()

    def get_bz_k_points(self):
        """Return all the k-points in the 1. Brillouin zone.

        The coordinates are relative to reciprocal latice vectors."""
        return self.driver.get_bz_k_points()

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        return self.driver.get_number_of_spins()

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        return self.driver.get_spin_polarized()

    def get_ibz_k_points(self):
        """Return k-points in the irreducible part of the Brillouin zone.

        The coordinates are relative to reciprocal latice vectors."""
        return self.driver.get_ibz_k_points()

    def get_k_point_weights(self):
        """Weights of the k-points.

        The sum of all weights is one."""
        return self.driver.get_k_point_weights()

    def get_pseudo_density(self, spin=None, pad=True, gather = False):
        """Return pseudo-density array.

        If *spin* is not given, then the total density is returned.
        Otherwise, the spin up or down density is returned (spin=0 or
        1)."""
        return self.driver.get_pseudo_density(spin=spin, pad=pad, gather=gather) / (units['Bohr'] ** 3)

    def get_effective_potential(self, spin=0, pad=True):
        """Return pseudo-effective-potential array."""
        return self.driver.get_effective_potential(spin=spin, pad=pad)

    def get_pseudo_wave_function(self, band=None, kpt=0, spin=0, broadcast=True,
                                 pad=True):
        """Return pseudo-wave-function array."""
        return self.driver.get_pseudo_wave_function(band=band, kpt=kpt, spin=spin, broadcast=broadcast, pad=pad)

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        return self.driver.get_eigenvalues(kpt=kpt, spin=spin)

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        return self.driver.get_occupation_numbers(kpt=kpt, spin=spin)

    def get_fermi_level(self):
        """Return the Fermi level."""
        return self.driver.get_fermi_level() * units['Ry']

    def initial_wannier(self, initialwannier, kpointgrid, fixedstates,
                        edf, spin, nbands):
        """Initial guess for the shape of wannier functions.

        Use initial guess for wannier orbitals to determine rotation
        matrices U and C.
        """
        raise NotImplementedError

    def get_wannier_localization_matrix(self, nbands, dirG, kpoint,
                                        nextkpoint, G_I, spin):
        """Calculate integrals for maximally localized Wannier functions."""
        raise NotImplementedError

    def get_magnetic_moment(self, atoms=None):
        """Return the total magnetic moment."""
        return self.driver.get_magnetic_moment(atoms=atoms)

    def get_number_of_grid_points(self, gather = True):
        """Return the shape of arrays."""
        return self.driver.get_number_of_grid_points(gather=gather)

def write2qe(outf, atoms, basefile = None, **kwargs):
    if basefile is None :
        ase.io.write(outf, atoms, format = 'espresso-in', **kwargs)
    else :
        ntyp = len(set(atoms.symbols))
        nat = len(atoms)
        with open(basefile, 'r') as fr :
            with open(outf, 'w') as fw :
                for line in fr :
                    if 'ntyp' in line :
                        x = line.index("=") + 1
                        line = line[:x] + ' ' + str(ntyp) + '\n'
                    elif 'nat' in line :
                        nat_old = int(line.split('=')[1])
                        x = line.index("=") + 1
                        line = line[:x] + ' ' + str(nat) + '\n'
                    elif 'cell_parameters' in line.lower() :
                        for i in range(3):
                            fr.readline()
                            line += '{0[0]:.14f} {0[1]:.14f} {0[2]:.14f}\n'.format(atoms.cell[i])
                    elif 'atomic_positions' in line.lower():
                        line = 'ATOMIC_POSITIONS angstrom\n'
                        for i in range(nat_old):
                            fr.readline()
                        for s, p in zip(atoms.symbols, atoms.positions):
                            line += '{0:4s} {1[0]:.14f} {1[1]:.14f} {1[2]:.14f}\n'.format(s, p)
                    fw.write(line)
