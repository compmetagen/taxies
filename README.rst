Example:

.. code-block:: sh

    TAXIES_DB_DIR=./db
    taxies genmark genome.fna $TAXIES_DB_DIR markers.fna -t 4
    taxies gentoph markers.fna $TAXIES_DB_DIR tophits.txt -t 4
    taxies gentax tophits.txt tax.txt

