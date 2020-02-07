Usage
=========


Pygmes is self documented in its help. Call it via:
.. code-block:: shell
    
    pygmes --help


To run pygmes on a fasta file please call:

.. code-block:: shell
    
    pygmes -i <input.fna> -o outdir --db database.dmnd  --ncores 16


Metagenomics mode
-----------------

To run pygmes in metagenomics mode you will need
all your bins in a single folder. They should all be 
from the same assembly, or of the same enviorement.

.. code-block:: shell
 
    pygmes -i <folder> -o outdir --db database.dmnd --meta --ncores 16

We recommend using 16 cores as this will speed up the analysis.
