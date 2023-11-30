#!/usr/bin/env python3
import argparse
from qepy.cui import qex

def get_args():
    parser = argparse.ArgumentParser(description='QEpy tasks :')
    #-----------------------------------------------------------------------
    for prog in qex.QEBINS :
        parser = qex.get_parse(prog=prog, parser=parser)
    #-----------------------------------------------------------------------
    args = parser.parse_args()
    return args

def main():
    import sys
    for prog in qex.QEBINS :
        key = '--' + prog
        if key in sys.argv :
            return qex.main(prog)
    else :
        get_args()
