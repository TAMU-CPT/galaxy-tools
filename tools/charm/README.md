# CHAdo Resource Manager

This is a re-write of charm to be a set of utilities for managing gmod product instances, rather than just chado. Mostly focused on a Tripal/Chado/JBrowse/WebApollo deployment

## Required Tools

- jbrowse
    - add new instances based on chado
    - update existing ones based on diff?
- bulkfiles
- chado
    - add-organism
    - gmod-bulk-load-gff3.pl
    - process.py
    - process-fasta.py
    - post-process-features.py
    - post-process-extract-chromosome.py
