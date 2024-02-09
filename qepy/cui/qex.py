#!/usr/bin/env python3
import argparse
from importlib import import_module
import qepy
from qepy.qebins import QEBINS
from qepy.io import set_input

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
    return args

def run(prog, qecmd=None, inputfile=None):
    if qecmd :
        qepy.qepy_modules.qepy_sys.set_command_line(qecmd)
    if inputfile is not None: set_input(inputfile)
    mod = import_module(QEBINS[prog][0])
    job = getattr(mod, QEBINS[prog][1])
    job()

def main(prog):
    args = get_args(prog)
    run(prog, qecmd=args.qecmd)
