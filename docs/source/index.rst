
pygmes
==================================


pygmes is a simple wrapper for GeneMark-ES. 

It will first run GeneMark-ES in fungal training mode on the used fasta. If this succeeds in will terminate there.

Else it will use a number of pre-trained models to predict proteins.
The protein file with the most amino acids will be compared againsta  protein database
using diamonds blastp to estimate a lineage. This lineage will be used
to find the best model from a number of 750 pre trained models.

In a final step the best taxonomic models (up to n) will be ran to predict proteins.
The protein prediction with the most amino acids in the prediction will
be used as the final output.

This way we can use GeneMark-ES on fragmented genomes that do not have enough data
to use the self training algorythm. Obviously this will impact the quality of 
the predicted proteins, but might be helpfull and sufficient for some
use cases.


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
    Usage <usage.rst>
    API <api.rst>




Indices and tables
==========================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
