.. _documentation:

=====================
Documentation
=====================

We use the reStructuredText_ markup language and the Sphinx_ to generate the documentation.

.. _Sphinx: http://www.sphinx-doc.org/en/master/
.. _reStructuredText: http://docutils.sourceforge.net/rst.html


Installing Sphinx and extensions
================================

.. highlight:: bash

Sphinx
------

Installation::

    $ pip install sphinx_rtd_theme --user

Extensions
----------

    * *nbsphinx* : Jupyter Notebook Tools for Sphinx::

        $ pip install nbsphinx --user

    * *sphinx-panels* : creating panels in a grid layout or as drop-downs::

        $ pip install sphinx-panels --user

Pandoc
------

`Installation (Ubuntu) <https://pandoc.org/installing.html>`_::

    $ sudo apt install pandoc
