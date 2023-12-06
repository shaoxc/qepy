import re
from importlib import import_module
from qepy import driver

optional=['qepy_cetddft']

def test_modules():
    file = driver.__file__
    pattern = re.compile(r'(qepy_.*?\..+?)\(')

    attrs = []
    with open(file, 'r') as fh:
        for line in fh:
            m = pattern.findall(line)
            attrs.extend(m)

    missing = []
    for item in attrs:
        l = item.split('.')
        try:
            mod = import_module('.'.join(l[:-1]))
        except Exception as e:
            if l[0] not in optional: raise e
            continue
        if not hasattr(mod, l[-1]):
            missing.append(item)

    if len(missing)>0:
        print('Following modules are not found:\n', '\n'.join(missing))
        raise AttributeError('Some modules are not found', len(missing))


if __name__ == "__main__":
    tests = [item for item in globals() if item.startswith('test_')]
    for func in sorted(tests):
        globals()[func]()
