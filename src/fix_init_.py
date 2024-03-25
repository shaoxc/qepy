import sys
from pathlib import Path

def fix_init(filename):
    init_new = """
# control the output
import sys
from importlib import import_module
from qepy.core import Logger, env
class QEpyLib :
    def __init__(self, **kwargs):
        qepylib = import_module(pname)
        sys.modules[pname] = self
        self.qepylib = qepylib

    def __getattr__(self, attr):
        attr_value = getattr(self.qepylib, attr)
        if '__array__' not in attr :
            attr_value = Logger.stdout2file(attr_value, fileobj=env['STDOUT'])
        return attr_value
qepylib = QEpyLib()
"""
    initfile = Path(filename)
    if initfile.is_file():
        mods = []
        with open(initfile, 'r+') as fh :
            lines = fh.readlines()
            fh.seek(0)
            if 'pname' not in lines[1]:
                pname = lines[1].split()[1]
                lines[1] = f"pname = '{pname}'\n" + init_new + lines[1]

                for line in lines :
                    if line.startswith('import qepy_'):
                        m = line.split()[-1]
                        v = m.partition('.')[2]
                        mods.append(v + ' = ' + m + '\n')
                    fh.write(line)
                fh.write('\n')
                for line in mods:
                    fh.write(line)
    return


if __name__ == "__main__":
    filename = '__init__.py'
    if len(sys.argv)>1:
        filename = sys.argv[1]
    fix_init(filename)
