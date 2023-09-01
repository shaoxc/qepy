import re
import os
from pathlib import Path

def replace_program_with_subroutine(contents):
    pattern = re.compile(r'program\s+(\w\w+)(.*?)end\s+program\s+\1', flags=re.IGNORECASE | re.DOTALL)
    pattern_only = re.compile(r'program\s+(\w\w+)', flags=re.IGNORECASE)
    has_contains = re.compile(r'\n\s*contains', flags=re.IGNORECASE)
    end_only = re.compile(r'(\n\s*end\s*\n)', flags=re.IGNORECASE)
    match = pattern.search(contents)
    if match is None :
        print(f"!WARN: Can not find end program in {filename}")
        match = pattern_only.search(contents)
        if match is None :
            raise ValueError(f"!!!ERROR: Can not find program in {filename}")
        else :
            results = pattern_only.sub(r'SUBROUTINE \1()', contents, count = 1)
            if not has_contains.search(results):
                m = end_only.search(results)
                if m :
                    results = end_only.sub(r'\nCONTAINS\n!\n', results, count = 1)
                    results += '\nEND SUBROUTINE'
    else :
        if has_contains.search(contents):
            results = pattern.sub(r'SUBROUTINE \1()\2END SUBROUTINE \1', contents, count = 1)
        else :
            # print('No contains', filename)
            results = pattern.sub(r'SUBROUTINE \1()\2CONTAINS\n!\n', contents, count = 1)
            results += f'\nEND SUBROUTINE {match.group(1)}'
    prog = match.group(1)
    return results, prog

def replace_value(contents, dicts):
    for k, v in dicts.items() :
        contents = re.sub(rf'\b{k}\b', rf'{v}', contents, flags=re.IGNORECASE, count=0)
    return contents

def comment_external(contents, value):
    pattern = None
    for v in value.split():
        pattern = re.compile(rf'\n(.*?external.*?{v})', flags=re.IGNORECASE)
        contents = pattern.sub(r'\n! \1', contents, count = 0)
    return contents

def replace_get_args(contents):
    sys_args = 'USE qepy_sys,           ONLY : COMMAND_ARGUMENT_COUNT, GET_COMMAND_ARGUMENT'
    pattern = re.compile(r'(\s*)(implicit none)', flags=re.IGNORECASE)
    contents = pattern.sub(rf'\1!!\1{sys_args}\1!!\1\2', contents, count = 0)
    return contents

def ini2files_qex(filename):
    from configparser import ConfigParser
    config = ConfigParser(allow_no_value=True)
    config.optionxform = str
    config.read(filename)
    #
    do_cmd = {Path(v).name :k for k, v in config['CMD'].items()}
    do_args = [Path(item).name for item in config['ARGS']]
    do_external = {Path(k).name :v for k, v in config['EXTERNAL'].items()}
    do_replace = {Path(k).name :v for k, v in config['REPLACE'].items()}
    qebins = {k:Path(v).name for k, v in config['CMD'].items()}
    files = set([v for k, v in config['CMD'].items()])
    files |= set(config['ARGS'])
    files |= set(config['EXTERNAL'].keys())
    files |= set(config['REPLACE'].keys())
    #
    path = os.environ.get('qedir', '')
    if path : path = path + os.sep

    for file in files :
        file = path + file
        name = Path(file).name
        with open(file, 'r') as fh:
            results = fh.read()
        if name in do_cmd :
            results, prog = replace_program_with_subroutine(results)
            for k, v in qebins.items():
                if v == name : qebins[k] = prog
        else :
            print(name)
        if name in do_args :
            results = replace_get_args(results)
        if name in do_external:
            results = comment_external(results, do_external[name])
        if name in do_replace :
            dicts = {}
            for item in do_replace[name].split():
                x = list(item.split(':'))
                dicts[x[0]] = x[1]
            results = replace_value(results, dicts)
        with open(name, 'w') as fh:
            fh.write(results)

    qebins['cell2ibrav.x'] = 'cell2ibrav'
    import pprint
    pprint.pprint(qebins)


if __name__ == '__main__':
    import sys
    filename = sys.argv[1] if len(sys.argv)>1 else 'qepy_qex.ini'
    ini2files_qex(filename)
