import numpy as np
from qepy.driver import Driver
from qepy.io import QEInput
from qepy import constants
from ase.calculators.calculator import Calculator
from ase.units import create_units
from ase.geometry import wrap_positions
import ase.io

units = create_units('2006').copy()
units['Bohr'] = constants.BOHR_RADIUS_SI * 1E10
units['Ry'] = constants.RYTOEV

class QEpyCalculator(Calculator):
    """QEpy Calculator for ASE.

    QEpy initialization depends on the input file. Here are some ways to generate input file :

         - Give a file `inputfile` or/and a dict `qe_options` will generate a temporary input file with new atomic information.
         - Give input file name `inputfile` and `from_file = True` will read the structures to atoms (ase.Atoms) from `inputfile` and initialize the QEpy with `inputfile`.
         - Give dict `ase_espresso` will generate a temporary input file with `ase.io.espresso.write_espresso_in <https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#ase.io.espresso.write_espresso_in>`__ (Only works for program PW).


    Parameters
    ----------
    atoms : ase.Atoms
        ASE atoms object
    inputfile : str
        Name of QE input file
    from_file : str
        Read the structure from inputfile, and use the inputfile for the QE Initialization.
    wrap : bool
        If True, wrap the atoms back to cell before update the ions.
    extrapolation : str or bool
        If False, every step will re-run the QE from file.
    ase_espresso : dict
        A dictionary with some parameters to generate PW input.
    qe_options: dict
        A dictionary with input parameters for QE to generate QE input file.
    comm : object
        Parallel communicator
    ldescf : bool
        If True, also print the scf correction term in the iteratively mode
    iterative : bool
        Iteratively run the scf or tddft
    task : str
        Task of the driver :

          - 'scf' : scf task
          - 'optical' : Optical absorption spectrum (TDDFT)
          - 'nscf' : read file from scf

    embed : object (Fortran)
        embed object for QE
    prefix : str
        prefix of QE input/output
    outdir : str
        The output directory of QE
    logfile : str, file or bool
        The screen output of QE. If it is a file descriptor, it should open with 'w+'.

          - None : Show on the screen.
          - str  : Save to the given file.
          - True : Save to the temporary file, see also :func:`qepy.driver.Driver.get_output`.


    .. note::

         - The bool `gather` parameter in functions means the data is gathered in root or distributed in all cpus.
         - if `inputfile` is not None. The ase_espresso should be set. Which contains :

             + input_data: dict
             + pseudopotentials: dict
             + kpts: (int, int, int) or dict
             + koffset: (int, int, int)
             + crystal_coordinates: bool

            More details, see `ase.io.espresso.write_espresso_in <https://wiki.fysik.dtu.dk/ase/ase/io/formatoptions.html#ase.io.espresso.write_espresso_in>`__.

    """

    implemented_properties=['energy', 'forces', 'stress']

    def __init__(self, atoms = None, inputfile = None, from_file = False, wrap = False, extrapolation = True,
            ase_espresso = None, qe_options = None, comm = None, ldescf = False, iterative = False,
            task = 'scf', embed = None, prefix = None, outdir = None, logfile = None, **kwargs):
        Calculator.__init__(self, atoms = atoms, **kwargs)
        self.optimizer = None
        self.atoms = atoms
        self.inputfile = inputfile
        self.wrap = wrap
        self.extrapolation = extrapolation
        self.from_file = from_file
        self.ase_espresso = ase_espresso
        self.qe_options = qe_options
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
            self.inputfile = 'input_tmp.in'
        #
        self.atoms_save = None
        self.driver = None
        self.iter = 0
        self.qeinput = QEInput()
        #
        self.qepy_options = {
                'comm' : comm,
                'ldescf' : ldescf,
                'iterative' : iterative,
                'task' : task,
                'embed' : embed,
                'logfile' : logfile,
                }
        self.qepy_options.update(kwargs)
        self.restart()

    def restart(self):
        self._energy = None
        self._forces = None
        self._stress = None

    def driver_initialise(self, atoms = None, prog = 'pw', **kwargs):
        atoms = atoms or self.atoms
        if self.wrap : atoms.wrap()
        if not self.from_file :
            if self.ase_espresso :
                ase.io.write(self.inputfile, atoms, format = 'espresso-in', **self.ase_espresso)
            else :
                self.qeinput.write_qe_input(self.inputfile, atoms=atoms, basefile=self.basefile,
                        qe_options=self.qe_options, prog=prog)
        self.driver = Driver(self.inputfile, **self.qepy_options)

    def update_atoms(self, atoms = None, **kwargs):
        atoms = atoms or self.atoms
        #
        if not self.extrapolation :
            self.stop()
            self.driver_initialise(atoms = atoms)
            return
        #
        if self.wrap :
            positions = wrap_positions(atoms.positions, atoms.cell)
        else :
            positions = atoms.positions

        positions = positions /units['Bohr']
        #
        lattice = atoms.cell /units['Bohr']
        lattice_old =self.driver.get_ions_lattice()
        if np.allclose(lattice, lattice_old): lattice = None
        #
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
        """See :func:`qepy.driver.Driver.stop`"""
        return self.driver.stop(what = 'all', **kwargs)

    @property
    def is_root(self):
        """See :func:`qepy.driver.Driver.is_root`"""
        return self.driver.is_root

    def get_density(self, gather = True):
        """See :func:`qepy.driver.Driver.get_density`"""
        return self.driver.get_density(gather=gather) / (units['Bohr'] ** 3)

    def get_wave_function(self, band=None, kpt=0):
        """See :func:`qepy.driver.Driver.get_wave_function`"""
        return self.driver.get_wave_function(band=band, kpt=kpt)

    def get_number_of_k_points(self):
        """See :func:`qepy.driver.Driver.get_number_of_k_points`"""
        return self.driver.get_number_of_k_points()

    def get_potential_energy(self, atoms = None, **kwargs):
        """See :func:`qepy.driver.Driver.get_potential_energy`"""
        self.update_optimizer(atoms)
        return self.driver.get_energy(**kwargs) * units['Ry']

    def get_forces(self, atoms = None, icalc = 0):
        """See :func:`qepy.driver.Driver.get_forces`"""
        if self.update_optimizer(atoms) or self._forces is None:
            self._forces = self.driver.get_forces(icalc=icalc)* units['Ry'] / units['Bohr']
        return self._forces

    def get_stress(self, atoms = None):
        """Return the stress tensor in the Voigt order (xx, yy, zz, yz, xz, xy)"""
        if self.update_optimizer(atoms) or self._stress is None:
            stress = self.driver.get_stress()
            self._stress = np.array(
                [stress[0, 0], stress[1, 1], stress[2, 2],
                 stress[1, 2], stress[0, 2], stress[0, 1]]) * -1 * units['Ry'] / (units['Bohr'] ** 3)
        return self._stress

    def get_number_of_bands(self):
        """See :func:`qepy.driver.Driver.get_number_of_bands`"""
        return self.driver.get_number_of_bands()

    def get_xc_functional(self):
        """See :func:`qepy.driver.Driver.get_xc_functional`"""
        return self.driver.get_xc_functional()

    def get_bz_k_points(self):
        """See :func:`qepy.driver.Driver.get_bz_k_points`"""
        return self.driver.get_bz_k_points()

    def get_number_of_spins(self):
        """See :func:`qepy.driver.Driver.get_number_of_spins`"""
        return self.driver.get_number_of_spins()

    def get_spin_polarized(self):
        """See :func:`qepy.driver.Driver.get_spin_polarized`"""
        return self.driver.get_spin_polarized()

    def get_ibz_k_points(self):
        """See :func:`qepy.driver.Driver.get_ibz_k_points`"""
        return self.driver.get_ibz_k_points()

    def get_k_point_weights(self):
        """See :func:`qepy.driver.Driver.get_k_point_weights`"""
        return self.driver.get_k_point_weights()

    def get_pseudo_density(self, spin=None, pad=True, gather = False):
        """See :func:`qepy.driver.Driver.get_pseudo_density`"""
        return self.driver.get_pseudo_density(spin=spin, pad=pad, gather=gather) / (units['Bohr'] ** 3)

    def get_effective_potential(self, spin=0, pad=True):
        """See :func:`qepy.driver.Driver.get_effective_potential`"""
        return self.driver.get_effective_potential(spin=spin, pad=pad)

    def get_pseudo_wave_function(self, band=None, kpt=0, spin=0, broadcast=True,
                                 pad=True):
        """See :func:`qepy.driver.Driver.get_pseudo_wave_function`"""
        return self.driver.get_pseudo_wave_function(band=band, kpt=kpt, spin=spin, broadcast=broadcast, pad=pad)

    def get_eigenvalues(self, kpt=0, spin=0):
        """See :func:`qepy.driver.Driver.get_eigenvalues`"""
        return self.driver.get_eigenvalues(kpt=kpt, spin=spin)

    def get_occupation_numbers(self, kpt=0, spin=0):
        """See :func:`qepy.driver.Driver.get_occupation_numbers`"""
        return self.driver.get_occupation_numbers(kpt=kpt, spin=spin)

    def get_fermi_level(self):
        """See :func:`qepy.driver.Driver.get_fermi_level`"""
        return self.driver.get_fermi_level() * units['Ry']

    # def initial_wannier(self, initialwannier, kpointgrid, fixedstates,
                        # edf, spin, nbands):
        # """See :func:`qepy.driver.Driver.initial_wannier`"""
        # raise NotImplementedError

    # def get_wannier_localization_matrix(self, nbands, dirG, kpoint,
                                        # nextkpoint, G_I, spin):
        # """See :func:`qepy.driver.Driver.get_wannier_localization_matrix`"""
        # raise NotImplementedError

    def get_magnetic_moment(self, atoms=None):
        """See :func:`qepy.driver.Driver.get_magnetic_moment`"""
        return self.driver.get_magnetic_moment(atoms=atoms)

    def get_number_of_grid_points(self, gather = True):
        """See :func:`qepy.driver.Driver.get_number_of_grid_points`"""
        return self.driver.get_number_of_grid_points(gather=gather)

    @staticmethod
    def get_atoms_from_qepy():
        """Return the atom.Atoms from QE."""
        atoms = Driver.get_ase_atoms()
        return atoms
