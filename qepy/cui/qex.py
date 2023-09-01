#!/usr/bin/env python3
import argparse
import qepy

QEBINS = {'alpha2f.x': 'alpha2f',
 'average.x': 'average',
 'bands.x': 'do_bands',
 'benchmark_libxc.x': 'benchmark_libxc',
 'cell2ibrav.x': 'cell2ibrav',
 'dist.x': 'pwscf',
 'dos.x': 'do_dos',
 'dynmat.x': 'dynmat',
 'epa.x': 'epa',
 'epsilon.x': 'epsilon',
 'ev.x': 'ev',
 'fd.x': 'fd',
 'fd_ef.x': 'fd_raman',
 'fd_ifc.x': 'fd_ifc',
 'fermi_proj.x': 'fermi_proj',
 'fermi_velocity.x': 'fermi_velocity',
 'fqha.x': 'fqha',
 'fs.x': 'fermisurface',
 'gww.x': 'gww',
 'gww_fit.x': 'gww_fit',
 'hp.x': 'hp_main',
 'ibrav2cell.x': 'ibrav2cell',
 'initial_state.x': 'initial_state',
 'kpoints.x': 'special_points',
 'lambda.x': 'elph',
 'ld1.x': 'ld1',
 'manypw.x': 'pwscf',
 'matdyn.x': 'matdyn',
 'molecularpdos.x': 'molecularpdos',
 'neb.x': 'neb',
 'open_grid.x': 'open_grid',
 'path_interpolation.x': 'images_interpolator',
 'ph.x': 'phonon',
 'plan_avg.x': 'plan_avg',
 'plotband.x': 'plotband',
 'plotproj.x': 'plotproj',
 'plotrho.x': 'plotrho',
 'pmw.x': 'pmw',
 'pp.x': 'pp',
 'ppacf.x': 'do_ppacf',
 'projwfc.x': 'do_projwfc',
 'pw.x': 'pwscf',
 'pw2bgw.x': 'pw2bgw',
 'pw2critic.x': 'pw2critic',
 'pw2gw.x': 'pw2gw',
 'pw2wannier90.x': 'pw2wannier90',
 'pw4gww.x': 'gwl_punch',
 'pwcond.x': 'pwcond',
 'q2qstar.x': 'Q2QSTAR',
 'q2r.x': 'q2r',
 'q2trans.x': 'q2trans',
 'q2trans_fd.x': 'q2trans_fd',
 'simple.x': 'simple',
 'simple_bse.x': 'simple_bse',
 'simple_ip.x': 'simple_ip',
 'sumpdos.x': 'sumpdos',
 'turbo_davidson.x': 'lr_dav_main',
 'turbo_eels.x': 'lr_eels_main',
 'turbo_lanczos.x': 'lr_main',
 'turbo_spectrum.x': 'lr_calculate_spectrum',
 'wannier_ham.x': 'wannier_ham',
 'wannier_plot.x': 'wannier_plot',
 'wfck2r.x': 'wfck2r'}

def get_parse(prog = 'pw.x', parser = None):
    if parser is None :
        parser = argparse.ArgumentParser(description=prog)
    parser.add_argument(f'--{prog}', dest=prog, action='store_true',
            default=False, help=f'Same as regular "{prog}" in QE.')
    return parser

def get_args(prog):
    parser = get_parse(prog)
    args, others = parser.parse_known_args()
    args.qecmd = ' '.join(others)
    qepy.qepy_sys.set_command_line(args.qecmd)
    return args

def run(args, prog):
    job = getattr(qepy, QEBINS[prog])
    print(job)
    job()

def main(prog):
    args = get_args(prog)
    run(args, prog)
