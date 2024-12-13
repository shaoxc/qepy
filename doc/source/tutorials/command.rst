.. _command:

==============
Run QEpy as QE
==============

QEpy can be used like a standard QE by providing input files and executing it directly from the command line.
Almost all QE packages are supported by this way.
To check the supported packages in QEpy:

.. code:: shell

    python -m qepy -h

Options starting with ``--`` and ending by ``.x`` correspond directly to the same package in the standard QE.
For example, you can run the ``pw.x`` package using:

.. code:: shell

    python -m qepy --pw.x -i qe.in

This is equivalent to running the regular QE command:

.. code:: shell

    pw.x -i qe.in

The parallel version QEpy also supports parallel execution with MPI. E.g.:

.. code:: shell

    mpirun -n 4 python -m qepy --pw.x -i qe.in

In addition to supporting standard QE packages, QEpy introduces new functionalities.
Options ending with ``.py`` are new implementations provided by QEpy.
For example, the ``--pw2pp.py`` package can analyze and convert the output of ``pw.x``.

To view the details of ``--pw2pp.py``, you can use:

.. code:: shell

    python -m qepy --pw2pp.py -h
