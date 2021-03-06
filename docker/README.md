# Taxies - Docker

Homepage: https://github.com/compmetagen/taxies.

This image includes:

 * Prodigal v2.6.3
 * HMMER v3.2.1
 * VSEARCH v2.9.1
 * Taxies (GitHub master)

## Quickstart

1. Download the latest version:

    `docker pull compmetagen/taxies`

2. Run an instance of the image, mounting the host working directory
    (e.g. ``/Users/davide/taxies``) on to the container working directory
    ``/taxies``:

    `docker run --rm -v /Users/davide/taxies:/taxies -w /taxies compmetagen/taxies taxies genmark ...`
    
    You need to write something like ``-v //c/Users/davide/taxies:/taxies`` if
    you are in Windows or ``-v /home/davide/taxies:/taxies`` in Linux. The
    ``--rm`` option automatically removes the container when it exits.
