.. highlight:: shell

==============
Advanced Usage
==============

Advanced blasting
_________________

Gene matching (*i.e.* blasting) between provided genomes and our internal database is performed with diamond_.

.. _diamond: https://github.com/bbuchfink/diamond

You can manually tweak the blasting options within diamond itself (at your own risk) as follows:

.. code-block:: console

    $ carve genome.faa --diamond-args="-e 1e-20 --top 20"

The default arguments are ``"--more-sensitive --top 10"``. Please see diamond's documentation for more details.


eggNOG-mapper
_____________


By default, **CarveMe** performs gene matching by *homology* search using diamond_.
However, you can also perform *orthology*-based search using eggNOG-mapper_.

.. _diamond: https://github.com/bbuchfink/diamond
.. _eggNOG-mapper: https://github.com/jhcepas/eggnog-mapper

For this you must first annotate your genome with eggNOG-mapper, and provide the output of eggNOG-mapper directly as
input to **CarveMe**:

.. code-block:: console

    $ carve --egg eggnog_output.tsv

Please make sure you install eggNOG-mapper from the *bigg* branch:

https://github.com/jhcepas/eggnog-mapper/tree/bigg


Media database
______________


**CarveMe** comes with a very small pre-built library_ of media compositions:

.. _library: https://github.com/cdanielmachado/carveme/blob/master/carveme/data/input/media_db.tsv

- ``LB`` (Lysogeny broth)
- ``LB[-O2]`` (Lysogeny broth, anaerobic)
- ``M9`` (Minimal M9 medium)
- ``M9[-O2]`` (Minimal M9 medium, anaerobic)
- ``M9[glyc]`` (Minimal M9 medium, glycerol as carbon source)

Additionally, you can provide your own media library for gap-filling:

.. code-block:: console

    $ carve genome.faa --gapfill X,Y,Z --mediadb mylibrary.tsv

The library must be a tab-separated file with four columns:

- *medium*: short id to be passed in command line (example: X)
- *description*: description of the medium (optional, example: Our magic X formula)
- *compound*: compound id (example: glc)
- *name*: compound name (optional, example: Glucose)

Please note that, at this moment, CarveMe only supports metabolite ids from the BiGG_ database.

.. _BiGG: http://bigg.ucsd.edu

Please feel free to contact us with suggestions of more media compositions to add to our default library.


SBML *flavor*
_____________

By default, **CarveMe** generates models compatible with the old **cobra toolbox** format.
This format is outdated but is still compatible with most constraint-based modeling tools.
The new format based on the **sbml-fbc2** specification is also supported.

You can specify your desired SBML *flavor* with the following flags:

.. code-block:: console

    $ carve genome.faa --cobra -o model.xml

    $ carve genome.faa --fbc2 -o model.xml


Ensemble modeling
_________________

Our model reconstruction algorithm is implemented as an MILP optimization problem. The generated model is structured
according to the solution to this problem. Often, one might want to explore how alternative solutions lead to slightly
different network structures, and consequently, predict different phenotypes.

**CarveMe** allows the generation of model ensembles. You only need to specify how many models you want to generate:

.. code-block:: console

    $ carve genome.faa -n 100 -o model.xml

This example would generate an ensemble of 100 models. Note that the ensemble is stored as a single SBML file, using
a compact notation (binary vectors) to represent the ensemble state of each reaction.

Some utility methods to read/write and perform simulations using ensemble models are implemented in *framed*.


Alternative universes
_____________________


**CarveMe** implements a top-down reconstruction approach that requires a well-curated universal model to be used as
template for the model *carving* process.

Currently, you can choose between the universal bacterial template, or two templates specialized for gram-positive and
gram-negative bacteria:

.. code-block:: console

    $ carve genome.faa -u grampos
    $ carve genome.faa -u gramneg

A script with some utility functions is available to help you build your own templates. For instructions please check:

.. code-block:: console

    $ build_universe -h

You can then provide your own customized universe model during reconstruction:

.. code-block:: console

    $ carve genome.faa --universe-file yeast_universe.xml


Experimental constraints
________________________

When you have experimental evidence for the presence/absence of a given set of reactions, you can provide this information
to improve the reconstruction process. According to the level of evidence, you can format your data as *soft* or *hard*
constraints. These can be applied to any kind of reaction present in the universe model (exchange, transport or enzymatic reactions).

**Soft constraints** are used to change the priority given to a set of reactions, as well as their expected direction.
They can be used when there is limited amount of evidence for some expected phenotype.
For instance, if the organism you are reconstructing is closely related to other organisms that are known to secrete a
given compound, you can include the respective exchange reaction as a *soft* constraint.

.. code-block:: console

    $ carve genome.faa --soft data.tsv

Where *data.tsv* is a tab-separated file with two columns, the reaction identifiers and the respective values.
Each value is one of the following: 1) reaction occurs in forward direction, -1) reaction occurs in backward direction,
0) reaction does not occur.

**Hard constraints** are used to force the fluxes through a given set of reactions during reconstruction. They can be used
when there is absolute evidence about a given phenotype. For instance, if you are reconstructing an obligatory anaerobe,
you can force the oxygen uptake rate to be zero.

.. code-block:: console

    $ carve genome.faa --hard data.tsv

Where *data.tsv* is a tab-separated file with three columns, reaction identifiers, lower bounds, and upper bounds.
Please use *hard* constraints with care, as they can make the reconstruction problem infeasible when incorrectly formulated.



