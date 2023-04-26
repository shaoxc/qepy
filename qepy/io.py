from collections import OrderedDict
import numpy as np

class QEInput(object):
    """Input for QE

    Read/write the QE input file. The input `qe_options` or the return of reading is a dictionary.
    There are two main types of QE input parameters : namelist and card.
    In the `qe_options`, for namelist, the key is always start with '&',  and the value is a dictionary.
    But for card, the key is the first line of the card, and the value is a list of string.

    *e.g.* ::

        qe_options = {
            '&control' : {
                'calculation' : "'scf'",
                'nstep' : 50,
                },
            'k_points automatic' : ['1 1 1 0 0 0'],
            }

    .. note ::

        If the type of value of namelist is string(character), the value also need contains with quotes.


    Parameters
    ----------
    filename : str
        Name of QE input file
    qe_options: dict
        A dictionary with input parameters for QE to generate QE input file.
    """

    def __init__(self, filename = None, qe_options = {}, **kwargs):
        self.filename = filename

        if self.filename is not None :
            self.qe_options = self.read_qe_input(self.filename, **kwargs)
            self.qe_options = self.update_options(options=qe_options, qe_options = self.qe_options)
        else :
            self.qe_options = qe_options

    def write_qe_input(self, filename, atoms = None, basefile = None, qe_options = {}, prog = 'pw', **kwargs):
        """Write the QE input file

        Parameters
        ----------
        filename : str
            The file name of output file
        atoms : ase.Atoms
            structure information with ase.Atoms
        basefile : str
            If set, will base on this file only update atomic information
        qe_options : dict
            Input parameters of QE
        """
        if basefile :
            options = qe_options
            qe_options = self.read_qe_input(basefile)
            qe_options = self.update_options(options=options, qe_options = qe_options, prog = prog)
        #
        if prog == 'pw' :
            if atoms is not None :
                qe_options = self.update_atoms(atoms, qe_options, prog=prog, **kwargs)
        if hasattr(filename, 'write'):
            fh = filename
        else :
            fh = open(filename, 'w')
        #
        options = QEOPTIONS.get(prog, {}).copy()
        options.update(qe_options)
        #
        for key, value in options.items() :
            if key.startswith('&'):
                fh.write(key.upper() + '\n')
                for k, v in value.items() :
                    fh.write('   ' + k + ' = ' + str(v) + '\n')
                fh.write('/\n\n')
            else :
                l = key.split()
                l[0] = l[0].upper()
                key = ' '.join(l)
                fh.write(key + '\n')
                for v in value :
                    fh.write(str(v) + '\n')
                fh.write('\n')

        if not hasattr(filename, 'write'): fh.close()

    @classmethod
    def update_atoms(cls, atoms, qe_options, prog = 'pw', **kwargs):
        """update atomic information

        Parameters
        ----------
        atoms : ase.Atoms
            atoms has all atomic information
        qe_options : dict
            Input parameters of QE
        """
        #
        if prog != 'pw' : return qe_options
        #
        if '&system' not in qe_options :
            qe_options['&system'] = {}
        qe_options['&system']['ibrav'] = 0
        qe_options['&system']['nat'] = len(atoms)
        ntyp = len(set(atoms.symbols))
        ntyp_old = len(qe_options.get('atomic_species', []))
        if ntyp > ntyp_old :
            raise ValueError("The number of 'ATOMIC_SPECIES' not fit the number of types of atoms.")
        qe_options['&system']['ntyp'] = ntyp_old
        # keys = [k.split()[0] for k in qe_options]
        keys = list(qe_options.keys())
        for key in keys :
            if key.startswith(('cell_parameters', 'atomic_positions')) :
                qe_options.pop(key, None)
        #
        key = 'cell_parameters angstrom'
        qe_options[key] = []
        for i in range(3):
            line = '{0[0]:.14f} {0[1]:.14f} {0[2]:.14f}'.format(atoms.cell[i])
            qe_options[key].append(line)
        #
        key = 'atomic_positions angstrom'
        qe_options[key] = []
        for s, p in zip(atoms.symbols, atoms.positions):
            line = '{0:4s} {1[0]:.14f} {1[1]:.14f} {1[2]:.14f}'.format(s, p)
            qe_options[key].append(line)
        #
        return qe_options

    @classmethod
    def update_options(cls, options = {}, qe_options = {}, **kwargs):
        """update options

        Use *options* to update *qe_options*

        Parameters
        ----------
        options : dict
            New options
        qe_options : dict
            QE options
        """

        options = options or {}
        qe_options = qe_options or {}
        keys = list(qe_options.keys())

        for k, v in options.items():
            k = k.lower()
            if k.startswith('&'):
                if k in qe_options :
                    qe_options[k].update(v)
                else :
                    qe_options[k] = v.copy()
            else :
                kc = k.split()[0]
                for item in keys :
                    if kc in item :
                        qe_options.pop(item)
                        break
                qe_options[k] = v
        # correct the qe_options
        for k, v in options.items():
            if 'cell_parameters' in k:
                qe_options['&system']['ibrav'] = 0
            elif 'atomic_positions' in k:
                qe_options['&system']['nat'] = len(v)
                ntyp = len(set([x.split()[0] for x in v]))
                ntyp_old = len(qe_options.get('atomic_species', []))
                if ntyp > ntyp_old :
                    raise ValueError("The number of 'ATOMIC_SPECIES' not fit the number of types of atoms.")
                qe_options['&system']['ntyp'] = ntyp_old
        return qe_options

    def read_qe_input(self, filename, **kwargs):
        """read the QE input file

        Parameters
        ----------
        filename : str
            file name of QE input file
        """
        options = {}
        with open(filename, 'r') as fr :
            fh = self.iter_lines(fr, sep='!', **kwargs)
            for line in fh :
                if line.startswith('&'):
                    options[line.split()[0].lower()] = self.read_namelist(fh, **kwargs)
                else :
                    fh.send(line)
                    break
            for line in fh :
                l = line.split()
                l[0] = l[0].lower()
                key = ' '.join(l)
                options[key] = self.read_card(fh, **kwargs)
        return options

    def read_namelist(self, fh, **kwargs):
        """read the namelist

        Parameters
        ----------
        fh : iter
            modified file handler
        """
        options = {}
        for line in fh :
            if line == '/' :
                break
            else :
                k, v = line.split('=')
                options[k.strip()] = v.strip(',').strip()
        else :
            raise ValueError("The namelist not closed with '/'.")
        return options

    def read_card(self, fh, **kwargs):
        """read the card

        The comments lines which start with '#' will be kept. Here assuming if all the characters of first word
        are alphabet letters and length is greater than 6, this line will be next card.

        Parameters
        ----------
        fh : iter
            modified file handler
        """
        l = []
        for line in fh :
            a = line.split()[0].replace('_', '')
            if a.isalpha() and len(a) > 6 : # assuming only card name has more than 6 characters
                fh.send(line)
                break
            else :
                l.append(line)
        return l

    def iter_lines(self, fh, sep='!', **kwargs):
        """return the non-comment part of non-empty line

        `sep` is the comment characters. The comment line or empty line will be skip.
        If the line is mixed, only the front non-comment part will return.

        Parameters
        ----------
        fh : object
            file handler
        sep : str, list or tuple
            The delimiter string, multiple delimiters can given by a list
        """
        if isinstance(sep, (tuple, list)):
            seps = sep[1:]
            sep = sep[0]
        else :
            seps = []
        for line in fh:
            for s in seps :
                line = line.replace(s,sep)
            line = line.split(sep)[0].strip()
            if len(line) > 0 :
                line = yield line
                if line is not None :
                    yield line
                    yield line


# All QE namelist from 7.1
QEOPTIONS={
    "all_currents" : OrderedDict.fromkeys(
        ["&energy_current"], {}),
    "bands" : OrderedDict.fromkeys(
        ["&bands"], {}),
    "bgw2pw" : OrderedDict.fromkeys(
        ["&input_bgw2pw"], {}),
    "cp" : OrderedDict.fromkeys(
        ["&control", "&system", "&electrons", "&ions", "&cell", "&press_ai", "&wannier", "atomic_species"], {}),
    "cppp" : OrderedDict.fromkeys(
        ["&inputpp"], {}),
    "davidson" : OrderedDict.fromkeys(
        ["&lr_input", "&lr_dav"], {}),
    "dos" : OrderedDict.fromkeys(
        ["&dos"], {}),
    "dynmat" : OrderedDict.fromkeys(
        ["&input"], {}),
    "eels" : OrderedDict.fromkeys(
        ["&lr_input", "&lr_control"], {}),
    "hp" : OrderedDict.fromkeys(
        ["&inputhp"], {}),
    "importexport_binary" : OrderedDict.fromkeys(
        ["&inputpp"], {}),
    "kcw" : OrderedDict.fromkeys(
        ["&control", "&wannier", "&screen", "&ham"], {}),
    "lanczos" : OrderedDict.fromkeys(
        ["&lr_input", "&lr_control", "&lr_post"], {}),
    "ld1" : OrderedDict.fromkeys(
        ["&input", "&inputp", "&test"], {}),
    "magnons" : OrderedDict.fromkeys(
        ["&lr_input", "&lr_control"], {}),
    "matdyn" : OrderedDict.fromkeys(
        ["&input"], {}),
    "molecularpdos" : OrderedDict.fromkeys(
        ["&inputmopdos"], {}),
    "neb" : OrderedDict.fromkeys(
        ["&path"], {}),
    "ph" : OrderedDict.fromkeys(
        ["&inputph"], {}),
    "postahc" : OrderedDict.fromkeys(
        ["&input"], {}),
    "pp" : OrderedDict.fromkeys(
        ["&inputpp", "&plot"], {}),
    "projwfc" : OrderedDict.fromkeys(
        ["&projwfc"], {}),
    "pw" : OrderedDict.fromkeys(
        ["&control", "&system", "&electrons", "&ions", "&cell", "&fcp", "&rism", "atomic_species"], {}),
    "pw2bgw" : OrderedDict.fromkeys(
        ["&input_pw2bgw"], {}),
    "pwcond" : OrderedDict.fromkeys(
        ["&inputcond"], {}),
    "q2r" : OrderedDict.fromkeys(
        ["&input"], {}),
    "spectrum" : OrderedDict.fromkeys(
        ["&lr_input"], {}),
    }


class QEOutput(object):
    def __init__(self, fh = None, nat = None, **kwargs):
        self.fh = fh
        self.nat = nat
        self.kwargs = kwargs

    @classmethod
    def get_forces_all(cls, fh = None, nat = None, **kwargs):
        if fh is None and hasattr(cls, 'fh') : fh = cls.fh
        if fh is None : raise ValueError('Please give a output "fh"')
        if nat is None and hasattr(cls, 'nat') : nat = cls.nat
        if nat is None : raise ValueError('Please give number of atoms "nat"')
        if isinstance(fh, list) : fh = iter(fh)
        forces={}
        lstart = False
        for line in fh:
            if len(line.strip())<3: continue
            if 'Total force' in line: break
            if 'Forces acting' in line:
                lstart = True
                line=next(fh)
                key = 'total'
            elif lstart:
                l = list(line.split())
                key=l[1] if l[0]=='The' else l[0]
            if lstart:
                fs = []
                for i in range(nat):
                    line=next(fh)
                    f=list(map(float, line.split()[-3:]))
                    fs.append(f)
                forces[key.lower()] = np.asarray(fs)
        return forces

    @classmethod
    def get_stress_all(cls, fh = None, **kwargs):
        if fh is None and hasattr(cls, 'fh') : fh = cls.fh
        if fh is None : raise ValueError('Please give a output "fh"')
        if isinstance(fh, list) : fh = iter(fh)
        stress ={}
        lstart = False
        for line in fh:
            if len(line.strip())<3: continue
            if lstart and 'stress' not in line: break
            if 'total   stress' in line:
                lstart = True
                line=next(fh)
                key = 'total'
                l = list(line.split())
            elif lstart:
                l = list(line.split())
                key = l[0]
            if lstart:
                ss = []
                for i in range(3):
                    s=list(map(float, line.split()[-3:]))
                    ss.append(s)
                    line=next(fh)
                stress[key.lower()] = np.asarray(ss)
        return stress

    @classmethod
    def get_forces(cls, fh = None, nat = None, **kwargs):
        return cls.get_forces_all(fh=fh, nat=nat, **kwargs)['total']

    @classmethod
    def get_stress(cls, fh = None, **kwargs):
        return cls.get_stress_all(fh=fh, **kwargs)['total']
