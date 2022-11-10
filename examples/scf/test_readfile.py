from qepy.driver import Driver

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

def main():
    driver = Driver(comm = comm, prefix = 'tmp', task = 'nscf')
    #
    energy = driver.get_energy()
    if driver.is_root :
        print('energy :\n', energy)
    driver.stop(what = 'no')


if __name__ == "__main__":
    main()
