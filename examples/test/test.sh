#!/usr/bin/env sh
check_exit() {
	if [ $1 -ne 0 ]
	then
		echo "!!!ERROR: something wrong in the test : exit code = $1"
		exit $1
	fi
}

serial='y'
parallel=''
if [ $# -gt 1 ]
then
	parallel='y'
elif [ $# -eq 1 ]
then
	if [ "$1" = "p" ]
	then
		serial=''
		parallel='y'
	fi
fi

if [ $serial ]; then
	echo "####################Test serial version####################"
	for f in *py
	do
		python3 -m pytest --cov-report term-missing --cov=./ $f
		check_exit $?
	done
fi

if [ $parallel ]; then
	echo "####################Test MPI version#######################"
	for f in *py
	do
		mpirun -n 2 python3 -m pytest --with-mpi $f
		check_exit $?
	done
fi
