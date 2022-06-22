from qepy.driver import Driver

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
except Exception:
    comm = None

def main():
    inputfile = './qe_in.in'
    driver = Driver(inputfile, comm, iterative = True)
    for i in range(60):
        driver.diagonalize()
        driver.mix()
        converged = driver.check_convergence()
        if converged : break
    #
    energy = driver.get_energy()
    forces = driver.get_forces()
    stress = driver.get_stress()
    if driver.is_root :
        print('converged :\n', converged)
        print('energy :\n', energy)
        print('forces :\n', forces)
        print('stress :\n', stress)
    driver.stop()


if __name__ == "__main__":
    main()
