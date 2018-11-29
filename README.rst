Example:

.. code-block:: sh

    TAXIES_DB_DIR=./db
    taxies genmark genome.fna $TAXIES_DB_DIR markers.fna -t 4
    taxies gentoph markers.fna $TAXIES_DB_DIR tophits.txt -t 4
    taxies gentax tophits.txt tax.txt


Docker image
^^^^^^^^^^^^
DockerHub: https://hub.docker.com/r/compmetagen/taxies/


Download:

.. code-block:: sh

    docker pull compmetagen/taxies

Build: 

.. code-block:: sh

    docker build docker/ -t taxies


Singularity image
^^^^^^^^^^^^^^^^^
SingularityHub: https://www.singularity-hub.org/collections/1963


Download:

.. code-block:: sh

    singularity pull --name taxies.simg shub://compmetagen/taxies

Build: 

.. code-block:: sh

    singularity build taxies.img Singularity