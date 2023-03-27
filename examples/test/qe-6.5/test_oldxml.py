import qepy
import numpy as np
import pathlib

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

path = pathlib.Path(__file__).resolve().parent / 'DATA/oldxml'
path = pathlib.Path(__file__).resolve().parent / 'DATA/oldxml_collect'

def test_oldxml():
    inputobj = qepy.qepy_common.input_base()
    inputobj.prefix = 'al_oldxml'
    inputobj.tmp_dir = str(path) + '/'
    if comm : inputobj.my_world_comm = comm.py2f()
    qepy.qepy_initial(inputobj)
    qepy.oldxml_read_file()
    embed = qepy.qepy_common.embed_base()
    qepy.qepy_calc_energies(embed)
    energy = embed.etotal
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    qepy.qepy_stop_run(0, what = 'no')

def test_oldxml_collect():
    inputobj = qepy.qepy_common.input_base()
    inputobj.prefix = 'al_oldxml'
    inputobj.tmp_dir = str(path) + '/'
    if comm : inputobj.my_world_comm = comm.py2f()
    qepy.qepy_initial(inputobj)
    qepy.oldxml_read_file()
    # 
    wfc = path / 'al_oldxml.wfc1'
    if wfc.is_file(): wfc.unlink()
    # 
    embed = qepy.qepy_common.embed_base()
    qepy.qepy_calc_energies(embed)
    energy = embed.etotal
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    qepy.qepy_stop_run(0, what = 'no')

if __name__ == "__main__":
    test_oldxml()
    test_oldxml_collect()
