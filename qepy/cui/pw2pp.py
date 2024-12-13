#!/usr/bin/env python3
import argparse
from qepy.driver import Driver
from qepy import qepy_pp

def get_parse(parser = None):
    if parser is None :
        parser = argparse.ArgumentParser(description='Analyze and convert from pw.x output')
    parser.add_argument('-i', '--input', dest='input', type=str, action='store',
            default=None, help='The input file of pw.x')
    parser.add_argument('-o', '--output', dest='output', type=str, action='store',
            default='qepy.pp', help='The output of file')
    parser.add_argument('--prefix', dest='prefix', type=str, action='store',
            default=None, help='prefix of files saved by program pw.x')
    parser.add_argument('--outdir', dest='outdir', type=str, action='store',
            default=None, help='directory of pw.x output')
    parser.add_argument('--plot_num', dest='plot_num', type=int, action='store',
            default=0, help='Selects what to save, default is electron density.')
    parser.add_argument('--calc_energy', dest='calc_energy', action='store_true',
            default=False, help='Calculate the energy from the output data.')
    return parser

def run():
    #
    parser = get_parse()
    args = parser.parse_args()
    inputfile = args.input
    output = args.output
    prefix = args.prefix
    outdir = args.outdir
    #
    if inputfile:
        driver = Driver(inputfile=inputfile)
        driver.pwscf_restart()
    elif prefix:
        driver = Driver(prefix = prefix, outdir=outdir, task = 'nscf')
    else:
        raise ValueError('Please provide a pw.x input file or "prefix" and "outdir"')

    if args.calc_energy:
        driver.calc_energy()

    filplot = output
    plot_num = args.plot_num
    ## These will go to the args
    kpoint = [0,0]
    kband = [0,0]
    spin_component = 0
    sample_bias = 0.01
    z = 1.0
    dz = 0.05
    lsign = 0
    emin = -999.0
    emax = +999.0
    ##
    qepy_pp.punch_plot(filplot, plot_num, sample_bias, z, dz,
            emin, emax, kpoint, kband, spin_component, lsign)

def main():
    run()


if __name__ == "__main__":
    main()
