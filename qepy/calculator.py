import numpy as np
import qepy
from ase.calculators.calculator import Calculator
from ase.units import create_units
import ase.io

units = create_units('2006')
units['Bohr'] = qepy.constants.BOHR_RADIUS_SI * 1E10
units['Ry'] = qepy.constants.RYTOEV

class QEpyCalculator(Calculator):
    """QEpy calculator for ase"""
    implemented_properties=['energy', 'forces', 'stress']

    def __init__(self, atoms = None, comm = None, task = 'scf', embed = None, inputfile = None, input_data = None, **kwargs):
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
            ase.io.write(self.inputfile, self.atoms, format = 'espresso-in', **self.parameters)
        if self.comm is None :
            comm = self.comm
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
            pos = atoms.positions.T / atoms.cell.cellpar()[0]
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
        qepy.qepy_mod.qepy_get_rho(self.density)
        return self.density / (units['Bohr'] ** 3)
