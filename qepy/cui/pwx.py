#!/usr/bin/env python3
import argparse
import qepy

def get_parse():
    parser = argparse.ArgumentParser(description='pw.x')
    parser.add_argument('--pw.x', dest='pw.x', action='store_true',
            default=False, help='Same as regular pw.x in QE.')
    return parser

def get_args():
    parser = get_parse()
    args, others = parser.parse_known_args()
    args.qecmd = 'qecmd ' + ' '.join(others)
    return args

def run(args):
    qepy.command_line_options.set_command_line_(args.qecmd)
    qepy.pwscf()

def main():
    args = get_args()
    run(args)
