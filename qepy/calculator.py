import numpy as np
import qepy
from ase.calculators.calculator import Calculator
from ase.units import create_units
from ase.geometry import wrap_positions
import ase.io

units = create_units('2006')
units['Bohr'] = qepy.constants.BOHR_RADIUS_SI * 1E10
units['Ry'] = qepy.constants.RYTOEV

class QEpyCalculator(Calculator):
    """QEpy calculator for ase"""
    implemented_properties=['energy', 'forces', 'stress']

    def __init__(self, atoms = None, comm = None, task = 'scf', embed = None, inputfile = None,
            input_data = None, wrap = False, **kwargs):
        Calculator.__init__(self, atoms = atoms, input_data = input_data, **kwargs)
        self.optimizer = None
        self.restart()
        self.task = task
        self.iter = 0
        if embed is None :
            embed = qepy.qepy_common.embed_base()
            self.olevel = 2
        else :
            self.olevel = 0
        self.embed = embed
        self.density = np.zeros((1, 1))
        self.atoms = atoms
        self.comm = comm
        self.inputfile = inputfile
        self.wrap = wrap
        if self.inputfile is not None :
            self.atoms = ase.io.read(self.inputfile, format='espresso-in')

    def restart(self):
        self._energy = None
        self._forces = None
        self._stress = np.zeros((3, 3), order='F')

    @property
    def rank(self):
        if self.comm is None :
            return 0
        else :
            return self.comm.rank

    def driver_initialise(self, **kwargs):
        if self.inputfile is None :
            self.inputfile = 'qepy_input.in'
            if self.wrap : self.atoms.wrap()
            ase.io.write(self.inputfile, self.atoms, format = 'espresso-in', **self.parameters)
        if self.comm is None :
            comm = None
        else :
            comm = self.comm.py2f()

        if self.task == 'optical' :
            qepy.qepy_tddft_main_initial(self.inputfile, comm)
            qepy.read_file()
        else :
            qepy.qepy_pwscf(self.inputfile, comm)

    def update_atoms(self, atoms = None, first = False, update = 0, **kwargs):
        atoms = atoms or self.atoms
        if not first :
            if self.wrap :
                positions = wrap_positions(atoms.positions, atoms.cell)
            else :
                positions = atoms.positions
            pos = positions.T / atoms.cell.cellpar()[0]
            qepy.qepy_mod.qepy_update_ions(self.embed, pos, update)
        if self.task == 'optical' :
            if first :
                self.embed.tddft.initial = True
                self.embed.tddft.finish = False
                self.embed.tddft.nstep = 900000 # Any large enough number
                qepy.qepy_tddft_readin(self.prefix + self._input_ext)
                qepy.qepy_tddft_main_setup(self.embed)

    def check_restart(self, atoms = None):
        restart = True

        if self.atoms :
            if atoms is None or atoms == self.atoms :
                restart = False

        if restart :
            if atoms is not None : self.atoms = atoms.copy()
            self.restart()
            if self.iter > 0 : self.update_atoms(self.atoms)

        if self.iter == 0 :
            restart = True
            self.driver_initialise()

        return restart

    def update_optimizer(self, atoms = None):
        if self.check_restart(atoms):
            self.iter += 1
            if self.task == 'optical' :
                qepy.qepy_molecule_optical_absorption(self.embed)
            else :
                qepy.qepy_electrons_scf(0, 0, self.embed)

    def get_potential_energy(self, atoms = None, **kwargs):
        self.update_optimizer(atoms)
        if self.olevel == 0 :
            qepy.qepy_calc_energies(self.embed)
        self._energy = self.embed.etotal
        return self._energy * units['Ry']

    def get_forces(self, atoms = None, icalc = 0):
        self.update_optimizer(atoms)
        qepy.qepy_forces(icalc)
        self._forces = qepy.force_mod.get_array_force().T
        return self._forces * units['Ry'] / units['Bohr']

    def get_stress(self, atoms = None):
        self.update_optimizer(atoms)
        qepy.stress(self._stress)
        stress = np.array(
            [self._stress[0, 0], self._stress[1, 1], self._stress[2, 2],
             self._stress[1, 2], self._stress[0, 2], self._stress[0, 1]])
        return stress * -1 * units['Ry'] / (units['Bohr'] ** 3)

    def get_density(self):
        nr = np.zeros(3, dtype = 'int32')
        qepy.qepy_mod.qepy_get_grid(nr)
        nspin = qepy.lsda_mod.get_nspin()
        if self.rank == 0 :
            if np.prod(nr) != self.density.shape[0] or nspin != self.density.shape[0] :
                self.density = np.empty((np.prod(nr), nspin), order = 'F')
        else :
            self.density = np.empty((1, nspin), order = 'F')
        qepy.qepy_mod.qepy_get_rho(self.density)
        return self.density / (units['Bohr'] ** 3)

    #ASE interface to DFT-calculators <--
    def get_number_of_bands(self):
        """Return the number of bands."""
        return qepy.wvfct.get_nbnd()

    def get_xc_functional(self):
        """Return the XC-functional identifier.

        'LDA', 'PBE', ..."""
        return qepy.funct.get_dft_short().decode("utf-8")

    def get_bz_k_points(self):
        """Return all the k-points in the 1. Brillouin zone.

        The coordinates are relative to reciprocal latice vectors."""
        return qepy.klist.get_array_xk()

    def get_number_of_spins(self):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        return qepy.lsda_mod.get_nspin()

    def get_spin_polarized(self):
        """Is it a spin-polarized calculation?"""
        return bool(qepy.lsda_mod.get_lsda())

    def get_ibz_k_points(self):
        """Return k-points in the irreducible part of the Brillouin zone.

        The coordinates are relative to reciprocal latice vectors."""
        return qepy.klist.get_array_xk()[:, :self.get_number_of_k_points()]

    def get_k_point_weights(self):
        """Weights of the k-points.

        The sum of all weights is one."""
        return qepy.klist.get_array_wk()[:self.get_number_of_k_points()]

    def get_pseudo_density(self, spin=None, pad=True, inone = False):
        """Return pseudo-density array.

        If *spin* is not given, then the total density is returned.
        Otherwise, the spin up or down density is returned (spin=0 or
        1)."""
        nr = self.get_number_of_grid_points(inone = inone)
        nspin = qepy.lsda_mod.get_nspin()
        density = np.empty((np.prod(nr), nspin), order = 'F')
        qepy.qepy_mod.qepy_get_rho(density, inone = inone)
        if spin is None :
            return density.sum(axis=1) / (units['Bohr'] ** 3)
        else :
            return density[:, spin] / (units['Bohr'] ** 3)

    def get_effective_potential(self, spin=0, pad=True):
        """Return pseudo-effective-potential array."""
        return qepy.scf.get_array_vrs()

    def get_pseudo_wave_function(self, band=None, kpt=0, spin=0, broadcast=True,
                                 pad=True):
        """Return pseudo-wave-function array."""
        qepy.qepy_mod.qepy_get_evc(kpt + 1)
        evc = qepy.wavefunctions.get_array_evc()
        if band is None :
            return evc
        else :
            return evc[:, band]

    def get_eigenvalues(self, kpt=0, spin=0):
        """Return eigenvalue array."""
        return qepy.wvfct.get_array_et()[:, kpt]

    def get_occupation_numbers(self, kpt=0, spin=0):
        """Return occupation number array."""
        return qepy.wvfct.get_array_wg()[:, kpt]

    def get_fermi_level(self):
        """Return the Fermi level."""
        return qepy.ener.get_ef()*units['Ry']

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
        return qepy.lsda_mod.get_magtot()

    def get_number_of_grid_points(self, inone = True):
        """Return the shape of arrays."""
        nr = np.zeros(3, dtype = 'int32')
        qepy.qepy_mod.qepy_get_grid(nr, inone)
        return nr
    #ASE interface to DFT-calculators -->

    def get_number_of_k_points(self):
        return qepy.klist.get_nks()
