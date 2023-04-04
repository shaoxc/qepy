import qepy
import numpy as np
import pathlib

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None


def test_oldxml():
    path = pathlib.Path(__file__).resolve().parent / 'DATA/oldxml'
    if comm and comm.size > 1 : return
    inputobj = qepy.qepy_common.input_base()
    inputobj.prefix = 'al_oldxml'
    inputobj.tmp_dir = str(path) + '/'
    if comm : inputobj.my_world_comm = comm.py2f()
    embed = qepy.qepy_common.embed_base()
    qepy.qepy_initial(inputobj, embed = embed)
    qepy.oldxml_read_file()
    qepy.qepy_calc_energies()
    energy = embed.etotal
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    qepy.qepy_stop_run(0, what = 'no')

def test_oldxml_collect():
    path = pathlib.Path(__file__).resolve().parent / 'DATA/oldxml_collect'
    inputobj = qepy.qepy_common.input_base()
    inputobj.prefix = 'al_oldxml'
    inputobj.tmp_dir = str(path) + '/'
    if comm : inputobj.my_world_comm = comm.py2f()
    embed = qepy.qepy_common.embed_base()
    qepy.qepy_common.set_embed(embed)
    qepy.qepy_initial(inputobj)
    qepy.oldxml_read_file()
    qepy.qepy_calc_energies()
    # 
    if qepy.io_global.get_ionode():
        for f in path.glob('al_oldxml.wfc*'):
            f.unlink()
    # 
    energy = embed.etotal
    assert np.isclose(energy, -552.93477389, atol = 1E-6)
    qepy.qepy_stop_run(0, what = 'no')

if __name__ == "__main__":
    test_oldxml()
    test_oldxml_collect()
