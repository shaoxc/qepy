#!/usr/bin/env python3
import argparse
import qepy
from qepy.qebins import QEBINS

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
