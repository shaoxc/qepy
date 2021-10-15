import numpy as np
import qepy

class QEpyDriver :
    def __init__(self, inputfile, comm = None, ldescf = False, **kwargs):
        qepy.qepy_pwscf(inputfile, comm)
        embed = qepy.qepy_common.embed_base()
        embed.ldescf = ldescf
        qepy.control_flags.set_niter(1)
        self.embed = embed
        self.iter = 0

    def diagonalize(self, print_level = 2, **kwargs):
        self.iter += 1
        self.embed.mix_coef = -1.0
        qepy.qepy_electrons_scf(print_level, 0, self.embed)

    def mix(self, mix_coef = 0.7, print_level = 2):
        self.embed.mix_coef = mix_coef
        qepy.qepy_electrons_scf(print_level, 0, self.embed)

    def check_convergence(self, **kwargs):
        return qepy.control_flags.get_conv_elec()

    def get_energy(self, **kwargs):
        return self.embed.etotal

    def get_forces(self, icalc = 0, **kwargs):
        qepy.qepy_forces(icalc)
        forces = qepy.force_mod.get_array_force().T
        return forces

    def get_stress(self, **kwargs):
        stress = np.ones((3, 3), order='F')
        qepy.stress(stress)
        return stress

    def stop(self, what = 'all', **kwargs):
        qepy.qepy_stop_run(0, what = what)
