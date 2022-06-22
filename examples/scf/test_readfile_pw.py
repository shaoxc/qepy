from qepy.driver import Driver

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

def main():
    inputfile = './qe_in.in'
    driver = Driver(inputfile, comm)
    driver.pwscf_restart()
    #
    energy = driver.get_energy()
    if driver.is_root :
        print('energy :\n', energy)
    driver.stop(what = 'no')


if __name__ == "__main__":
    main()
