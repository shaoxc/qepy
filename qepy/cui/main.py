#!/usr/bin/env python3
import argparse
from importlib import import_module
from qepy.cui import qex

QEPYBINS = {
    'pw2pp.py' : 'qepy.cui.pw2pp',
}

def get_parse():
    def formatter(prog):
        return argparse.HelpFormatter(prog, max_help_position=30, width=100)
    parser = argparse.ArgumentParser(description='QEpy tasks :', formatter_class=formatter)
    #-----------------------------------------------------------------------
    for prog in qex.QEBINS :
        p = qex.get_parse(prog=prog)
        parser.add_argument(f'--{prog}', dest=prog, action='store_true',
                default=False, help = p.description)
    for prog, modfile in QEPYBINS.items():
        mod = import_module(modfile)
        p = mod.get_parse()
        parser.add_argument(f'--{prog}', dest=prog, action='store_true',
                default=False, help = p.description)
    #-----------------------------------------------------------------------
    return parser

def get_job_name(argv):
    engine = key = None
    for arg in argv[1:] :
        if arg.startswith('--') :
            key = arg[2:]
            if key in qex.QEBINS:
                engine = 'qe'
                break
            elif key in QEPYBINS:
                engine = 'qepy'
                break
    if engine: argv.remove('--'+key)
    return engine, key


def main():
    import sys
    engine, prog = get_job_name(sys.argv)
    if engine == 'qe' :
        return qex.main(prog)
    elif engine == 'qepy' :
        mod = import_module(QEPYBINS[prog])
        return mod.main()
    else:
        parser = get_parse()
        parser.print_help()
