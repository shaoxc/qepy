import numpy as np
import tempfile
import os
from functools import wraps
import qepy
from qepy.core import env
from qepy.io import QEInput
from qepy import constants

def gathered(function):
    @wraps(function)
    def wrapper(self, **kwargs):
        out = kwargs.get('out', None)
        gather = kwargs.get('gather', True)
        if out is None : out = self.create_array(gather=gather, kind='rho')
        if gather and self.nproc > 1 :
            arr = self.create_array(gather=False, kind='rho')
        else :
            arr = out
        #
        kwargs['out'] = arr
        results = function(self, **kwargs)
        #
        if gather and self.nproc > 1 :
            qepy.qepy_mod.qepy_get_value(arr, out, gather = True)
            if isinstance(results, (tuple, list)):
                results = out, *results[1:]
            else :
                results = out
        return results
    return wrapper

class Driver(object) :
    """
    The driver of QEpy.

    Parameters
    ----------
    inputfile : str
        Name of QE input file
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

    outdir : str
        The output directory of QE
    qe_options : dict
        A dictionary with input parameters for QE to generate QE input file.
    prog : str
        The name of QE program, default is `pw` which is pw.x in QE.
    progress : bool
        If True, will continue run the QE without clean the workspace, most times used for TDDFT after scf.
    kwargs : dict
        Other options


    .. note::

         - The bool `gather` parameter in functions means the data is gathered in root or distributed in all cpus.
         - Also contains some interface of ASE calculator, but the units are different.

             + positions : Bohr
             + lattice : Bohr
             + energy : Ry
             + forces : Ry/Bohr
             + stress : Ry/Bohr**3
             + potential : Ry
             + density : 1.0/Bohr**3

    """

    def __init__(self, inputfile = None, comm = None, ldescf = False, iterative = False,
             task = 'scf', embed = None, prefix = None, outdir = None, logfile = None,
             qe_options = None, prog = 'pw', progress = False, atoms = None, **kwargs):
        if embed is None :
            embed = qepy.qepy_common.embed_base()
        self.task = task
        self.embed = embed
        self.comm = comm
        self.inputfile = inputfile
        self.iterative = iterative
        self.prefix = prefix
        self.outdir = outdir
        self.logfile = logfile
        self.qe_options = qe_options
        self.prog = prog
        self.progress = progress
        self.atoms = atoms
        #
        self.embed.ldescf = ldescf
        #
        self.comm = comm
        self.qepy = qepy
        self.qeinput = QEInput()
        #
        self.driver_initialize()

    def _init_log(self):
        """_initialize the QE output."""
        self.fileobj_interact = False
        if self.logfile is not None :
            # qepy.qepy_mod.qepy_set_stdout(self.logfile)
            if isinstance(self.logfile, bool) and self.logfile :
                self.fileobj_interact = True
                self.fileobj = tempfile.NamedTemporaryFile('w+')
            elif hasattr(self.logfile, 'write'):
                self.fileobj = self.logfile
            else :
                self.fileobj = open(self.logfile, 'w+')
        else :
            self.fileobj = None
        env['STDOUT'] = self.fileobj
        return self.fileobj

    @property
    def comm(self):
        return self._comm

    @comm.setter
    def comm(self, value):
        self._comm = value
        if self.comm is not None and hasattr(self.comm, 'py2f') :
            self.commf = self.comm.py2f()
        else :
            self.commf = self.comm

    @property
    def is_root(self):
        """Whether it is the root.

        If the comm is set, root is rank == 0.
        Otherwise root is ionode of QE.
        """
        if hasattr(self.comm, 'rank'):
            return self.comm.rank == 0
        else :
            return qepy.io_global.get_ionode()

    @property
    def nproc(self):
        if hasattr(self.comm, 'size'):
            return self.comm.size
        else :
            return qepy.mp_world.get_nproc()

    @property
    def qe_is_mpi(self):
        return qepy.qepy_common.get_is_mpi()

    @property
    def qe_is_openmp(self):
        return qepy.qepy_common.get_is_openmp()

    def restart(self, prog=None, **kwargs):
        prog = prog or self.prog
        self.driver_initialize()

    def driver_initialize(self, **kwargs):
        """ Initialize the driver

        Parameters
        ----------
        inputfile : str
            Name of QE input file
        commf : object
            Parallel communicator (Fortran)
        iterative : bool
            Iteratively run the scf or tddft
        task : str
            Task of the driver :

              - 'scf' : scf task
              - 'optical' : Optical absorption spectrum (TDDFT)
              - 'nscf' : read file from scf

        qe_options: dict
            A dictionary with input parameters for QE to generate QE input file.

        """
        # stop the last driver and save new driver
        if not self.progress :
            if hasattr(env['DRIVER'], 'stop'): env['DRIVER'].stop()
            env['DRIVER'] = self
        #
        self._init_log()
        qepy.qepy_common.set_embed(self.embed)
        #
        inputfile=self.inputfile
        commf=self.commf
        task=self.task
        qe_options = self.qe_options
        prog = self.prog
        #
        if qe_options :
            inputfile, basefile = 'input_tmp.in', inputfile
            self.qeinput.write_qe_input(inputfile, atoms = self.atoms, basefile=basefile, qe_options=qe_options, prog=prog)
        #
        if task == 'optical' :
            self.tddft_initialize(inputfile=inputfile, commf = commf, **kwargs)
        elif task == 'nscf' :
            inputobj = qepy.qepy_common.input_base()
            if self.prefix : inputobj.prefix = self.prefix
            if self.outdir : inputobj.tmp_dir = str(self.outdir) + '/'
            if commf : inputobj.my_world_comm = commf
            qepy.qepy_initial(inputobj)
            tmpdir = inputobj.tmp_dir.decode().strip() + inputobj.prefix.decode().strip() + '.save' + '/'
            if os.path.isfile(tmpdir + 'data-file.xml'): # only works for qe-6.5 !
                if not hasattr(qepy, 'oldxml_read_file') :
                    raise AttributeError("Please reinstall the QEpy with 'oldxml=yes'.")
                qepy.oldxml_read_file()
            else :
                qepy.read_file()
            qepy.qepy_mod.qepy_open_files()
        else :
            qepy.qepy_pwscf(inputfile, commf)
            self.embed.iterative = self.iterative
            if self.embed.iterative :
                qepy.control_flags.set_niter(1)
        #
        self.density = np.zeros((1, 1))
        self.iter = 0

    def tddft_initialize(self, inputfile = None, commf = None, **kwargs):
        """ Initialize the tddft

        Parameters
        ----------
        inputfile : str
            Name of QE input file, which also contains `&inputtddft` section.
        commf : object
            Parallel communicator (Fortran)
        """
        if inputfile is None : inputfile = self.inputfile
        if commf is None : commf = self.commf
        #
        if self.progress :
            qepy.wvfct.get_array_g2kin()
            qepy.qepy_tddft_readin(inputfile)
        else :
            qepy.qepy_tddft_main_initial(inputfile, commf)
            qepy.read_file()
        qepy.qepy_tddft_main_setup()
        self.embed.tddft.iterative = self.iterative

    def diagonalize(self, print_level = 2, **kwargs):
        """Diagonalize the hamiltonian

        Parameters
        ----------
        print_level :
            The level of output of QE
        """
        self.iter += 1
        if self.task == 'optical' :
            qepy.qepy_molecule_optical_absorption()
        else :
            self.embed.mix_coef = -1.0
            qepy.qepy_electrons_scf(print_level, 0)

    def mix(self, mix_coef = 0.7, print_level = 2):
        """Mixing the density

        Parameters
        ----------
        print_level :
            The level of output of QE
        """
        if self.task == 'optical' : return
        self.embed.mix_coef = mix_coef
        qepy.qepy_electrons_scf(print_level, 0)

    def check_convergence(self, **kwargs):
        """Check the convergence of the SCF"""
        converged = bool(qepy.control_flags.get_conv_elec())
        if converged and not self.embed.initial : self.end_scf()
        return converged

    def get_scf_error(self, **kwargs):
        """Return the error of the scf"""
        if self.embed.iterative :
            return self.embed.dnorm
        else :
            return qepy.control_flags.get_scf_error()

    def get_scf_steps(self, **kwargs):
        """Return the number of steps of scf"""
        if self.embed.iterative :
            return self.iter
        else :
            return qepy.control_flags.get_n_scf_steps()

    def scf(self, print_level = 2, maxiter = None, original = False, **kwargs):
        """Run the scf/tddft until converged or maximum number of iterations"""
        if maxiter is not None and not self.embed.iterative :
            qepy.control_flags.set_niter(maxiter)
        if self.task == 'optical' :
            qepy.qepy_molecule_optical_absorption()
        elif not self.embed.iterative and self.embed.exttype < 2 :
            # Use electrons to support hybrid xc functional
            return self.electrons(original=original)
        else :
            qepy.qepy_electrons_scf(print_level, 0)
        return self.embed.etotal

    def non_scf(self, **kwargs):
        # fix some saved variables from last scf calculations
        qepy.control_flags.set_lscf(0)
        qepy.control_flags.set_lbfgs(0)
        qepy.control_flags.set_lmd(0)
        qepy.control_flags.set_lwf(0)
        #
        qepy.non_scf()
        return qepy.ener.get_etot()

    def electrons(self, original = False, **kwargs):
        if original :
            qepy.electrons()
        else :
            qepy.qepy_electrons()
        return qepy.ener.get_etot()

    def end_scf(self, **kwargs):
        """End the scf and clean the scf workspace. Only need run it in iterative mode"""
        if self.embed.iterative :
            self.embed.finish = True
            qepy.qepy_electrons_scf(0, 0)

    def stop(self, exit_status = 0, what = 'all', print_flag = 0, **kwargs):
        """stop.

        Parameters
        ----------
        exit_status : int
            0 : QE will remove it temporary files
        print_flag : int
            0 : Not print time informations
        what : str
             see :func:`qepy.driver.Driver.save`.
        """
        if self.task == 'optical' :
            self.tddft_stop(exit_status, print_flag = print_flag, what = what, **kwargs)
        else :
            qepy.qepy_stop_run(exit_status, print_flag = print_flag, what = what, finalize = False)

        if hasattr(self.fileobj, 'close'): self.fileobj.close()
        qepy.qepy_clean_saved()
        #
        env['DRIVER'] = None
        env['STDOUT'] = None
        #

    def tddft_restart(self, istep=None, **kwargs):
        """Restart the tddft from previous interrupted run.

        Parameters
        ----------
        istep : int
            Start number of steps, just for output.
        """
        qepy.qepy_tddft_mod.qepy_cetddft_wfc2rho()
        if istep is not None :
            self.embed.tddft.istep = istep

    def tddft_stop(self, exit_status = 0, what = 'no', print_flag = 0, **kwargs):
        if self.embed.tddft.iterative :
            self.embed.tddft.finish = True
            qepy.qepy_molecule_optical_absorption()
        #! Do not save the PW files, otherwise the initial wfcs will be overwritten.
        qepy.qepy_stop_run(exit_status, print_flag = print_flag, what = 'no', finalize = False)
        qepy.qepy_stop_tddft(exit_status)

    def save(self, what = 'all', **kwargs):
        """
        Save the QE data to the disk

        Parameters
        ----------
        what : str

          * what = 'all' : write xml data file, charge density, wavefunctions
                           (for final data);
          * what = 'config' : write xml data file and charge density; also,
                              for nks=1, wavefunctions in plain binary format
                              (see why in comments below). For intermediate
                              or incomplete results;
          * what = 'config-nowf' : write xml data file and charge density only;
                           (save density);
          * what = 'config-init' : write xml data file only excluding final results
                                   (for dry run, can be called at early stages).

          see PW/src/punch.f90
        """
        qepy.punch(what)

    def update_run_options(self, qe_options = {}, **kwargs):
        pass

    def get_energy(self, **kwargs):
        """Return the total energy."""
        if abs(qepy.ener.get_etot()) > 1E-16 :
            energy = qepy.ener.get_etot()
        elif abs(self.embed.etotal) > 1E-16 :
            energy = self.embed.etotal
        else :
            energy = self.calc_energy(**kwargs)
        return energy

    def calc_energy(self, **kwargs):
        """Calculate the energy with the pw2casino of QE."""
        qepy.qepy_calc_energies()
        return self.embed.etotal

    def update_ions(self, positions = None, lattice = None, update = 0, **kwargs):
        """update the ions of QE

        Parameters
        ----------
        positions : np.ndarray (n, 3)
            positions of ions
        lattice : 3x3 matrix
            lattice of ions
        update : int
             - 0 : update everything (default)
             - 1 : only update atomic information, not update potential

        """
        positions = positions.T
        if lattice is not None :
            lattice = lattice.T
            if not qepy.cellmd.get_lmovecell():
                raise ValueError(" Lattice update only works for variable-cell simulations.\n Please restart the QEpy with calculation= 'vc-relax' or 'vc-md'")
            qepy.qepy_mod.qepy_update_ions(positions, update, lattice)
        else :
            qepy.qepy_mod.qepy_update_ions(positions, update)

    def pwscf_restart(self, oldxml=False, starting_pot='file', starting_wfc='file'):
        """Read PW ouput/restart files.

        Parameters
        ----------
        oldxml : bool
             - True : read the old format xml file (qe<6.4)
             - False : read the new format xml file (qe>6.3)
        """
        if oldxml :
            if not hasattr(qepy, 'oldxml_pw_restart') :
                raise AttributeError("Please reinstall the QEpy with 'oldxml=yes'.")
            qepy.oldxml_pw_restart.pw_readfile('header')
            qepy.oldxml_pw_restart.pw_readfile('reset')
            qepy.oldxml_pw_restart.pw_readfile('dim')
            qepy.oldxml_pw_restart.pw_readfile('bs')
            if qepy.basis.get_starting_pot().strip() != starting_pot :
                qepy.basis.set_starting_pot(starting_pot)
                qepy.oldxml_potinit()
            if qepy.basis.get_starting_wfc().strip() != starting_wfc :
                qepy.basis.set_starting_wfc(starting_wfc)
                qepy.oldxml_wfcinit()
        else :
            qepy.qepy_mod.qepy_restart_from_xml()
            if qepy.basis.get_starting_pot().strip() != starting_pot :
                qepy.basis.set_starting_pot(starting_pot)
                qepy.potinit()
            if qepy.basis.get_starting_wfc().strip() != starting_wfc :
                qepy.basis.set_starting_wfc(starting_wfc)
                qepy.wfcinit()

    def create_array(self, gather = True, kind = 'rho'):
        """Return an empty array in real space."""
        if kind == 'rho' :
            nspin = qepy.lsda_mod.get_nspin()
            if gather and self.nproc > 1 :
                nr = self.get_number_of_grid_points(gather = gather)
                if self.is_root :
                    out = np.zeros((np.prod(nr), nspin), order = 'F')
                else :
                    out = np.zeros((1, nspin), order = 'F')
            else :
                nnr = self.embed.dfftp.nnr
                out = np.zeros((nnr, nspin), order = 'F')
        else :
            raise ValueError('Only support "rho"')
        return out

    def get_density(self, gather = True, out = None):
        """Return density array in real space."""
        if out is None : out = self.create_array(gather=gather, kind='rho')
        qepy.qepy_mod.qepy_get_rho(out, gather = gather)
        return out

    def get_kinetic_energy_density(self, gather = True, out = None):
        """Return density array in real space."""
        if out is None : out = self.create_array(gather=gather, kind='rho')
        qepy.qepy_mod.qepy_get_tau(out, gather = gather)
        return out

    def get_wave_function(self, band=None, kpt=0):
        """Return wave-function array in real space."""
        qepy.qepy_mod.qepy_get_evc(kpt + 1)
        nrs = np.zeros(3, dtype = 'int32')
        qepy.qepy_mod.qepy_get_grid_smooth(nrs)
        if self.is_root :
            wf = np.empty(np.prod(nrs), order = 'F', dtype = np.complex128)
        else :
            wf = np.empty(1, order = 'F', dtype = np.complex128)
        if band is None :
            band = np.arange(self.get_number_of_bands())
        else :
            band = np.asarray(band)
            if band.ndim == 0 : band = [band]
        wfs = []
        for ibnd in band :
            qepy.qepy_mod.qepy_get_wf(kpt + 1, ibnd + 1, wf)
            wfs.append(wf.copy())
        return wfs

    def get_dipole_tddft(self):
        """Return the total dipole of tddft task."""
        # dipole = qepy.qepy_tddft_common.get_array_dipole().copy()
        dipole = self.embed.tddft.dipole
        return dipole

    def set_external_potential(self, potential, exttype = None, gather = True, extene  = None, **kwargs):
        """Set the external potential.

        Parameters
        ----------
        potential : (nnr, nspin)
            The external potential
        exttype : int
            The type of external potential

                 ==== ============================== ===
                 type potential                      bin
                 ==== ============================== ===
                  0   external                       000
                  1   pseudo                         001
                  2   hartree                        010
                  3   pseudo + hartree               011
                  4   xc                             100
                  5   pseudo + xc                    101
                  6   hartree + xc                   110
                  7   pseudo + hartree + xc          111
                 ==== ============================== ===

        """
        if exttype is not None :
            self.embed.exttype = exttype
        if extene is not None :
            self.embed.extene = extene
        else :
            self.embed.extene = 0.0
        #
        if potential is None : potential = np.zeros((1, 1))
        if potential.ndim != 2 : raise ValueError("The array should be 2-d.")
        #
        qepy.qepy_mod.qepy_set_extpot(potential, gather = gather)

    def get_output(self):
        """Return the output of QE.

        It depends on the `logfile` of the initialization of the driver :

          - None : return None.
          - str  : return all outputs.
          - True : It will return the output from last time.
        """
        if self.fileobj is not None :
            if self.fileobj_interact :
                self.fileobj.flush()
                self.fileobj.seek(0)
                lines = self.fileobj.readlines()
                self.fileobj.close()
                #
                self._init_log()
                #
                return lines
            else :
                self.fileobj.seek(0)
                return self.fileobj.readlines()
        else :
            return None

    @gathered
    def get_elf(self, gather = True, out = None, **kwargs):
        """Return electron localization function."""
        qepy.do_elf(out)
        return out

    @gathered
    def get_rdg(self, gather = True, out = None, **kwargs):
        """Return electron localization function."""
        qepy.do_rdg(out)
        return out
#-----------------------------------------------------------------------

    @gathered
    def get_local_pp(self, gather = True, out = None, **kwargs):
        """Return local component of the pseudopotential."""
        for i in range(out.shape[1]):
            out[:, i] = qepy.scf.get_array_vltot()
        return out

    @gathered
    def get_hartree(self, gather = True, out = None, add=False, **kwargs):
        """Return hartree information."""
        if not add : out[:] = 0.0
        ehart, charge = qepy.v_h(self.embed.rho.of_g[:,0], out)
        return out, ehart, charge

    @gathered
    def get_exchange_correlation(self, gather = True, out = None, tau=None, **kwargs):
        """Return exchange-correlation information.

        TODO :
            The interface will changed in the new version of QE!

        Note :
            For metaGGA, only return scattered tau
        """
        etxc = vtxc = 0.0
        rho_obj = self.embed.rho
        rho_core = qepy.scf.get_array_rho_core()
        rhog_core = qepy.scf.get_array_rhog_core()
        if qepy.funct.dft_is_meta():
            if tau is None : tau = out*0.0
            qepy.v_xc_meta(rho_obj, rho_core, rhog_core, etxc, vtxc, out, tau)
            return out, etxc, vtxc, tau
        else :
            etxc, vtxc = qepy.v_xc(rho_obj, rho_core, rhog_core, out)
            return out, etxc, vtxc

    @gathered
    def get_density_functional(self, gather = True, out = None, add = False, **kwargs):
        """Return effective potential information.

        Note :
            Then final potential is saved in the v_obj
        """
        rho_obj = self.embed.rho
        rho_core = qepy.scf.get_array_rho_core()
        rhog_core = qepy.scf.get_array_rhog_core()
        v_obj = self.embed.v
        etotefield = 0.0
        #
        v_obj.of_r[:] = 0.0
        ehart, etxc, vtxc, eth, charge = qepy.qepy_v_of_rho(rho_obj, rho_core, rhog_core, etotefield, v_obj)
        info = [ehart, etxc, vtxc, eth, etotefield, charge]
        if add :
            out += v_obj.of_r
        else :
            out[:] = v_obj.of_r
        return out, *info

    def get_hartree_potential(self, gather = True, out = None, **kwargs):
        return self.get_hartree(gather=gather, out=out, **kwargs)[0]

    def get_exchange_correlation_potential(self, gather = True, out = None, **kwargs):
        return self.get_exchange_correlation(gather=gather, out=out, **kwargs)[0]

    def get_density_functional_potential(self, gather = True, out = None, **kwargs):
        return self.get_density_functional(gather=gather, out=out, **kwargs)[0]

    def get_effective_potential(self, gather = True, out = None, **kwargs):
        out = self.get_local_pp(gather=gather, out=out, **kwargs)
        out = self.get_density_functional_potential(gather=gather, out=out, add = True, **kwargs)
        return out
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
    #ASE Calculator
    def get_potential_energy(self, **kwargs):
        """Return the potential energy."""
        return self.get_energy()

    def get_forces(self, icalc = 0, **kwargs):
        """Return the total forces. (n, 3)

        Parameters
        ----------
        icalc : int

            ===== ============================= ===
            icalc without                       bin
            ===== ============================= ===
              0   all                           000
              1   no ewald                      001
              2   no local                      010
              3   no ewald + local              011
              4   no nlcc                       100
              5   no ewald + nlcc               101
              6   no local + nlcc               110
              7   no ewald + local + nlcc       111
            ===== ============================= ===
        """
        qepy.qepy_forces(icalc)
        forces = qepy.force_mod.get_array_force().T
        return forces

    @classmethod
    def get_stress(cls, **kwargs):
        """Return the stress (3, 3)."""
        stress = np.zeros((3, 3), order='F')
        qepy.stress(stress)
        return stress

    #ASE DFTCalculator
    @classmethod
    def get_number_of_bands(cls):
        """Return the number of bands."""
        return qepy.wvfct.get_nbnd()

    @classmethod
    def get_xc_functional(cls):
        """Return the XC-functional identifier.

        'LDA', 'PBE', ..."""
        return qepy.funct.get_dft_short().decode("utf-8")

    @classmethod
    def get_bz_k_points(cls):
        """Return all the k-points in the 1. Brillouin zone.

        The coordinates are relative to reciprocal latice vectors."""
        return qepy.klist.get_array_xk()

    @classmethod
    def get_number_of_spins(cls):
        """Return the number of spins in the calculation.

        Spin-paired calculations: 1, spin-polarized calculation: 2."""
        return qepy.lsda_mod.get_nspin()

    @classmethod
    def get_spin_polarized(cls):
        """Is it a spin-polarized calculation?"""
        return bool(qepy.lsda_mod.get_lsda())

    @classmethod
    def get_ibz_k_points(cls):
        """Return k-points in the irreducible part of the Brillouin zone.

        The coordinates are relative to reciprocal latice vectors."""
        xk = qepy.klist.get_array_xk()[:, :cls.get_number_of_k_points()].T
        return xk @ qepy.cell_base.get_array_at()

    @classmethod
    def get_k_point_weights(cls):
        """Weights of the k-points.

        The sum of all weights is one."""
        return qepy.klist.get_array_wk()[:cls.get_number_of_k_points()]

    def get_pseudo_density(self, spin=None, pad=True, gather = True):
        """Return pseudo-density array.

        If *spin* is not given, then the total density is returned.
        Otherwise, the spin up or down density is returned (spin=0 or
        1)."""
        density = self.get_density(gather = gather)
        if spin is None :
            return density.sum(axis=1)
        else :
            return density[:, spin]

    @classmethod
    def get_pseudo_wave_function(cls, band=None, kpt=0, spin=0, broadcast=True,
                                 pad=True):
        """Return pseudo-wave-function array."""
        qepy.qepy_mod.qepy_get_evc(kpt + 1)
        evc = qepy.wavefunctions.get_array_evc()
        if band is None :
            return evc
        else :
            return evc[:, band]

    @classmethod
    def get_eigenvalues(cls, kpt=0, spin=0):
        """Return eigenvalue array."""
        return qepy.wvfct.get_array_et()[:, kpt]

    @classmethod
    def get_occupation_numbers(cls, kpt=0, spin=0):
        """Return occupation number array."""
        return qepy.wvfct.get_array_wg()[:, kpt]

    @classmethod
    def get_fermi_level(cls):
        """Return the Fermi level."""
        return qepy.ener.get_ef()

    # def initial_wannier(self, initialwannier, kpointgrid, fixedstates,
                        # edf, spin, nbands):
        # """Initial guess for the shape of wannier functions.

        # Use initial guess for wannier orbitals to determine rotation
        # matrices U and C.
        # """
        # raise NotImplementedError

    # def get_wannier_localization_matrix(self, nbands, dirG, kpoint,
                                        # nextkpoint, G_I, spin):
        # """Calculate integrals for maximally localized Wannier functions."""
        # raise NotImplementedError

    @classmethod
    def get_magnetic_moment(cls, **kwargs):
        """Return the total magnetic moment."""
        return qepy.lsda_mod.get_magtot()

    @classmethod
    def get_number_of_grid_points(cls, gather = True):
        """Return the shape of arrays."""
        nr = np.zeros(3, dtype = 'int32')
        qepy.qepy_mod.qepy_get_grid(nr, gather)
        return nr
    #ASE DFTCalculator END
#-----------------------------------------------------------------------

    @classmethod
    def get_number_of_k_points(cls):
        """Return the number of kpoints."""
        return qepy.klist.get_nks()

    @classmethod
    def get_volume(cls):
        """Return the volume."""
        return qepy.cell_base.get_omega()

    @classmethod
    def get_ecutrho(cls):
        """Return the cutoff for density."""
        return qepy.gvect.get_ecutrho()

    @classmethod
    def get_ions_lattice(cls):
        """Return the lattice of ions."""
        alat = qepy.cell_base.get_alat()
        lattice = qepy.cell_base.get_array_at() * alat
        return lattice.T

    @classmethod
    def get_ions_positions(cls):
        """Return the cartesian positions of ions."""
        alat = qepy.cell_base.get_alat()
        pos = qepy.ions_base.get_array_tau().T * alat
        return pos

    @classmethod
    def get_ions_symbols(cls):
        """Return the symbols of ions."""
        ityp = qepy.ions_base.get_array_ityp() - 1
        nat = qepy.ions_base.get_nat()
        ntyp = qepy.ions_base.get_nsp()
        label = qepy.ions_base.get_array_atm().T.view('S3')[:,0].astype('U3')[:ntyp]
        label = [x.strip() for x in label]
        symbols = []
        for i in range(nat):
            symbols.append(label[ityp[i]])
        return symbols

    @classmethod
    def data2field(cls, data, cell = None, grid = None, rank = None):
        """QE data to dftpy DirectField.
        If data is None or small temporary array will return np.zeros(1).
        """
        from dftpy.field import DirectField
        from dftpy.grid import DirectGrid
        #
        if data is None or data.size < 8 : return np.zeros(1)
        #
        if cell is None : cell = cls.get_ions_lattice()
        if grid is None : grid = DirectGrid(lattice=cell, nr=cls.get_number_of_grid_points(), ecut=cls.get_ecutrho())
        if not rank :
            rank = data.shape[1] if data.ndim == 2 else 1
        #
        if data.size != grid.nnrR*rank :
            raise ValueError('The size of data not match the grid,\
                    please check the data or set another grid', data.size, grid.nnrR*rank)
        #
        field = DirectField(grid=grid, data=data.ravel(order='C'), order='F', rank=rank)
        return field

    @classmethod
    def field2data(cls, field, data = None):
        """dftpy DirectField to QE data.
        If the input field is not 3d array, will return np.zeros((1, 1)).
        """
        #
        if field is None or field.ndim < 3 : return np.zeros((1, 1))
        #
        ns = field.shape[0] if field.ndim == 4 else 1
        if data is None :
            nnrR = np.prod(cls.get_number_of_grid_points())
            data = np.empty((nnrR, ns))
        #
        if data.size != field.size :
            raise ValueError('The size of field is different with data,\
                    please check the field or set another data', field.size, data.size)
        #
        if ns == 1 : field = [field]
        for i in range(ns):
            data[:, i] = field[i].ravel(order = 'F')
        return data

    @classmethod
    def get_dftpy_grid(cls, nr = None, cell = None, mp = None, **kwargs):
        """Return the dftpy DirectGrid from QE."""
        from dftpy.grid import DirectGrid
        if cell is None : cell = cls.get_ions_lattice()
        if nr is None :
            nr = cls.get_number_of_grid_points()
            ecut = cls.get_ecutrho()
        else :
            ecut = None
        grid = DirectGrid(lattice=cell, nr=nr, ecut=ecut, mp=mp, **kwargs)
        return grid

    @classmethod
    def get_ase_atoms(cls):
        """Return the atom.Atoms from QE."""
        from ase.atoms import Atoms
        units_Bohr = constants.BOHR_RADIUS_SI * 1E10
        #
        symbols = cls.get_ions_symbols()
        positions = cls.get_ions_positions() * units_Bohr
        lattice = cls.get_ions_lattice() * units_Bohr
        atoms = Atoms(symbols = symbols, positions = positions, cell = lattice)
        return atoms

    @classmethod
    def get_dftpy_ions(cls):
        """Return the dftpy.Ions from QE."""
        from dftpy.ions import Ions
        atoms = cls.get_ase_atoms()
        return Ions.from_ase(atoms)

    @classmethod
    def set_density(cls, density, gather = True, **kwargs):
        """Set density array in real space."""
        #
        if density is None : density = np.zeros((1, 1))
        if density.ndim != 2 : raise ValueError("The array should be 2-d.")
        #
        qepy.qepy_mod.qepy_set_rho(density, gather = gather)

    @classmethod
    def switch_nlpp(cls, nhm=0, nbetam=0, nkb=0, nh=None, **kwargs):
        nhm_ = qepy.uspp_param.get_nhm()
        nbetam_ = qepy.uspp_param.get_nbetam()
        nh_ = qepy.uspp_param.get_array_nh().copy()
        nkb_ = qepy.uspp.get_nkb()

        if nh is None: nh = nh_ * 0

        qepy.uspp_param.set_nhm(nhm)
        qepy.uspp_param.set_nbetam(nbetam)
        qepy.uspp.set_nkb(nkb)
        qepy.uspp_param.set_array_nh(nh)

        pp_options = {
            'nhm' : nhm_,
            'nbetam' : nbetam_,
            'nh': nh_,
            'nkb': nkb_,
        }

        return pp_options

    @classmethod
    def update_exchange_correlation(cls, xc=None, libxc=None, **kwargs):
        if libxc : xc = None
        qepy.qepy_mod.qepy_set_dft(xc)
