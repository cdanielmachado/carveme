.. highlight:: shell

============
Installation
============

CarveMe can be easily installed using the **pip** package manager:

.. code-block:: console

    $ pip install carveme

Additionally, you must manually install the diamond_ sequence aligner:

.. code-block:: console

    $ conda install -c bioconda diamond



Solver
------

You also need an optimization solver. The recommended option is to install either IBM CPLEX_ or Gurobi_.
You need a full version with the respective academic license (demo version will not work). 

After installing CPLEX, **do not forget** to install the CPLEX python API.

**Alternatively**, you can use the open source SCIP_ solver. This will be slower (expect at least 10 min per execution).

.. code-block:: console

    $ conda install --c conda-forge pyscipopt

Also, remember to select the correct solver during execution:

.. code-block:: console

    $ carve [ARGS] --solver scip

.. _diamond: https://github.com/bbuchfink/diamond
.. _SCIP: https://scipopt.org
.. _CPLEX: https://www.ibm.com/analytics/cplex-optimizer
.. _Gurobi: https://www.gurobi.com/downloads/gurobi-software
.. _documentation: https://www.ibm.com/support/knowledgecenter/SSSA5P_12.7.1/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html

Everything should be ready now! See the next section for instructions on how to start *carving*.
