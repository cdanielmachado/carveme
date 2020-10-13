.. highlight:: shell

=====
Usage
=====

Building a model
----------------

**CarveMe** provides a very simple command line interface to build models.
The most basic utilization is:

.. code-block:: console

    $ carve genome.faa

This will build a genome-scale metabolic model from the **genome** file.

By default **CarveMe** expects a protein FASTA file. Alternatively, you can also provide DNA sequences instead:

.. code-block:: console

    $ carve --dna genome.fna

Note that raw genome files are not supported. The FASTA file must be divided into individual genes.

It is possible to specify a different name or directory for the output file:

.. code-block:: console

    $ carve genome.faa --output model.xml

Short version:

.. code-block:: console

    $ carve genome.faa -o model.xml

If you want to produce a compressed SBML file, just change the extension (this is automatically supported by libSBML):

.. code-block:: console

    $ carve genome.faa -o model.xml.gz

Rather than providing the genome data yourself, you can also provide an NCBI RefSeq accession code.
This will automatically download the sequence and build the model:

.. code-block:: console

    $ carve --refseq GCF_000005845.2 -o ecoli_k12_mg1655.xml

If you have downloaded multiple genome sequences, you can run *recursive mode* to build multiple models in one call.

This will launch multiple parallel processes, which can decrease the overall computation time if you are running
**CarveMe** in a multi-core CPU or in a computing cluster:

.. code-block:: console

    $ carve -r myfolder/*.faa

This can be combined with *-o* to change the output folder:

.. code-block:: console

    $ carve -r myfolder/*.faa -o mymodels/


Gap Filling
-----------

**CarveMe** tries to predict the uptake and secretion capabilities of an organism only from genetic evidence,
and will produce a simulation-ready model without gap-filling for any particular media.

However, there are situations where you want to guarantee that the model is able to reproduce growth in one, or several,
experimentally verified media.

For instance, you can ensure the model reproduces growth on M9 and LB media:

.. code-block:: console

    $ carve genome.faa --gapfill M9,LB

Short version:

.. code-block:: console

    $ carve genome.faa -g M9,LB

Please see the *Advanced Usage* section on how to provide your own media compositions.

If you already have a model, and you just want to gap-fill it, you can do it with the *gapfill* utility function:

.. code-block:: console

    $ gapfill model.xml -m M9 -o new_model.xml

Please note that the result is not the same if you gap-fill during reconstruction. When you gap-fill during
reconstruction, the gene annotation scores are used to prioritize the reactions selected for gap-filling based on
genetic evidence. If you invoke *gapfill* alone, all potential gap-filling reactions are treated equally.

Finally, it is important to note that the models generated with **CarveMe** are not initialized with any
medium composition.

You can define the growth environment of the organism for simulation purposes by setting the flux bounds
of the exchange reactions yourself to match the respective medium composition.

Alternatively, you can tell **CarveMe** you want the model to come with a pre-defined medium composition.

.. code-block:: console

    $ carve genome.faa --init M9

Short version:

.. code-block:: console

    $ carve genome.faa -i M9

Note that this will not gap-fill the model, but only define the external environment for simulation purposes.

To simultaneously gap-fill and initialize the model for a desired medium, you must combine both flags:

.. code-block:: console

    $ carve genome.faa -g M9 -i M9

You are now a basic user. Happy *carving*!


Microbial Communities
---------------------

**CarveMe** enables the generation of microbial community models from single species models.

The most basic usage is:

.. code-block:: console

    $ merge_community organism_1.xml organism_2.xml ... organism_N.xml -o community.xml

or more simply:

.. code-block:: console

    $ merge_community *.xml -o community.xml

This generates an SBML file with a community where each organism is assigned to its own compartment and
a common community biomass equation is also generated. You can import the merged model into any simulation tool, just
as any normal constraint-based model and apply different types of simulation methods (FBA, FVA, etc...).
You can initialize the community with a pre-defined medium (just like during single-species reconstruction):

.. code-block:: console

    $ merge_community [input files] -i M9


