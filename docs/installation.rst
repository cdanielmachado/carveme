.. highlight:: shell

============
Installation
============

**CarveMe** requires Python 2.7, which is available for all major operating systems. We recommend the `Anaconda python
distribution <https://www.continuum.io/downloads>`_.

It can be easily installed using the **pip** package manager:

.. code-block:: console

    $ pip install carveme

This will automatically install other dependencies as well:

- framed_
- pandas_

.. _framed: https://github.com/cdanielmachado/framed
.. _pandas: https://pandas.pydata.org/

Additionally, you must manually install two external dependencies:

- diamond_
- IBM CPLEX_ Optimizer

.. _diamond: https://github.com/bbuchfink/diamond
.. _CPLEX: https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/

Note that you will need to register with IBM to obtain an academic license for CPLEX.

Finally, you should run the initialization script. This will download larger files
to build an internal database for *CarveMe*, which does not come bundled with the installer:

.. code-block:: console

    $ carveme_init

Everything should be ready now! See the next section for instructions on how to start *carving*.
