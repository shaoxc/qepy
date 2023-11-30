#!/usr/bin/env python
# coding: utf-8

import sys
from pathlib import Path
from f90wrap import parser as fparse
from f90wrap import fortran as ft
import multiprocessing
from multiprocessing import Pool
from functools import partial

def get_procedures(i, files=[], nprocs=1):
    n = len(files)//nprocs+1
    tree = fparse.read_files(files[i*n:i*n+n])
    return get_proces(tree)

def get_proces(tree, return_all=False):
    mods = []
    subs = []
    proces = []
    progs_files = []
    for node in ft.walk(tree):
        if isinstance(node, ft.Program):
            progs_files.append(node.filename)
        elif not isinstance(node, ft.Procedure):
            continue
        elif node.filename in progs_files: # skip the program files
            continue
        else:
            mod = ft.find_procedure_module(tree, node)
            subs.append(node)
            mods.append(mod)
            if mod is None:
                proces.append(node)
            elif mod not in proces:
                proces.append(mod)
    if return_all:
        return proces, progs_files, subs, mods
    else:
        return proces

def get_duplicate(proces):
    names = [item.name for item in proces]
    index_dict = {}
    duplicates = {}
    for i, item in enumerate(names):
        if item in index_dict:
            ip = index_dict[item]
            if proces[ip].filename == proces[i].filename: continue
            if item in duplicates:
                if proces[duplicates[item][1]].filename == proces[i].filename: continue
                duplicates[item].append(i)
            else:
                duplicates[item] = [ip, i]
        else:
            index_dict[item] = i
    return duplicates


path = sys.argv[1]
skips = ['external/', 'tests/', 'test.f90']

src = Path(path)
lists = src.glob('**/*.f90')
files = []
for item in lists:
    item = str(item)
    for s in skips:
        if s in item: break
    else:
        files.append(item)

nprocs = multiprocessing.cpu_count()
job = partial(get_procedures, files=files, nprocs=nprocs)
pool = Pool(processes=nprocs)
proces_all = pool.map(job, range(nprocs))
proces = sum(proces_all, [])

duplicates = get_duplicate(proces)

for k, v in duplicates.items():
    print(f'Duplicate functions: {k:<}')
    for i in v:
        print(f'\t {proces[i].filename}')
