# CHAdo Resource Manager

This is a re-write of charm to be a set of utilities for managing gmod product instances, rather than just chado. Mostly focused on a Tripal/Chado/JBrowse/WebApollo deployment

## Required Tools

- jbrowse
    - [ ] add new instances based on chado
    - [ ] update existing ones based on diff?
- bulkfiles
    - [ ] bulkfiles wrapper
- chado
    - [x] add-organism
    - [x] gmod-bulk-load-gff3.pl
    - [x] process.py
    - [x] process-fasta.py
    - [x] post-process-features.py
    - [x] post-process-extract-chromosome.py
