.. module:: qepy.calculator

The QEpyCalculator object
=========================


The :class:`QEpyCalculator` object is a calculator for ASE_::

  from qepy.calculator import QEpyCalculator
  driver = QEpyCalculator(inputfile='qe_in.in')

.. note::
   The units of `ASE Calculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html>`__ is different from :class:`qepy.driver.Driver` object.

.. _ASE: https://gitlab.com/ase/ase

List of all Methods
-------------------

.. autoclass:: QEpyCalculator
   :members:
