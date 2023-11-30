import re
from importlib import import_module
from qepy import driver


def test_modules():
    file = driver.__file__
    pattern = re.compile(r'(qepy\..+?)\(')

    attrs = []
    with open(file, 'r') as fh:
        for line in fh:
            m = pattern.findall(line)
            attrs.extend(m)

    missing = []
    for item in attrs:
        l = item.split('.')
        mod = import_module('.'.join(l[:-1]))
        if not hasattr(mod, l[-1]):
            missing.append(item)

    if len(missing)>0:
        print('Following modules are not found:\n', '\n'.join(missing))
        raise AttributeError('Some modules are not found', len(missing))


if __name__ == "__main__":
    test_modules()
