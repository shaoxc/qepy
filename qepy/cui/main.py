#!/usr/bin/env python3
import argparse

commands = {
        '--pw.x' : 'qepy.cui.pwx',
        }

def get_args():
    parser = argparse.ArgumentParser(description='pw.x')
    #-----------------------------------------------------------------------
    parser.description = 'QEpy tasks :\n' + '\n\t'.join(commands.keys())
    parser.add_argument('--pw.x', dest='pw.x', action='store_true',
            default=False, help='Same as regular pw.x in QE.')
    #-----------------------------------------------------------------------
    args = parser.parse_args()
    return args


def main():
    import sys
    from importlib import import_module
    for job in commands :
        if job in sys.argv :
            module = import_module(commands[job])
            command = getattr(module, 'main')
            return command()
    else :
        get_args()
