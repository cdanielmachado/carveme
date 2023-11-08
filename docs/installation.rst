.. highlight:: shell

============
Installation
============

CarveMe can be easily installed using the **pip** package manager:

.. code-block:: console

    $ pip install carveme

Additionally, you must manually install these dependencies:

- diamond_ (conda install -c bioconda diamond)
- SCIP_ solver [optional] (conda install --c conda-forge pyscipopt)
- IBM CPLEX_ or Gurobi_ optimizer [recommended] (full version with academic license)

.. _diamond: https://github.com/bbuchfink/diamond
.. _SCIP: https://scipopt.org
.. _CPLEX: https://www.ibm.com/analytics/cplex-optimizer
.. _Gurobi: https://www.gurobi.com/downloads/gurobi-software

Note that you will need to register with IBM to obtain an `academic license for CPLEX <https://www.ibm.com/academic/home>`_.

**IMPORTANT:** After installing CPLEX, do not forget to install the CPLEX python API (see the CPLEX documentation_ for details).

.. _documentation: https://www.ibm.com/support/knowledgecenter/SSSA5P_12.7.1/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html

Everything should be ready now! See the next section for instructions on how to start *carving*.
