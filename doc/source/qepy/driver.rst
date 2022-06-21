.. module:: qepy.driver

The Driver object
=================


The :class:`Driver` object is a simple way to call the QE/QEpy from python side. Here is how to define a simple scf job::

  from qepy import Driver
  driver = Driver('qe_in.in')

Here, the first argument specifies the name of `QE input file <https://www.quantum-espresso.org/Doc/INPUT_PW.html>`__.

List of all Methods
-------------------

.. autoclass:: Driver
   :members:
