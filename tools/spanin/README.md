# Overview

## The goal of this tool is to take putative protein sequences and determine if any of them are quality candidates for a potential spanin pair (OSP/ISP)
### Current Scripts:
    * `generate-putative-isp.py`
        * INPUT : Genomic FASTA
        * OUTPUT : Putative candidates for isp
    * `generate-putative-osp.py`
        * INPUT : Genomic FASTA
        * OUTPUT : Putative candidates for osp
    * `spaninFuncs.py`
        * Functions that drive the scripts within this directory

## Requirements:
    * regex
    * cpt.py

## Script Descriptions / Methodologies:
    * `generate-putative-isp.py`
        * INPUT : Genomic FASTA
        * METHODOLOGY 
            * Uses the `OrfFinder` object and locates __ALL__ potential start sequences, based on `ttg` / `atg` / `gtg`
            * This generates a set of output files consisting of FASTA / BED / GFF3 
            * With the aaFASTA, we read in each potential sequence and determine if it has a lipobox based off of: 
                * Min distance from start codon to lipobox
                * Max distance from start codon to lipobox
            * Within that range, a regular expression (RegEx) is used to determine if a lipobox is found.
        * OUTPUT : All matching description/sequence pairs are stored, and returned as a `putative_osp.fa` file.
    * `generate-putative-osp.py`
        * INPUT : Genomic FASTA
        * METHODOLOGY
        * Uses the `OrfFinder` object and locates __ALL__ potential start sequences, based on `ttg` / `atg` / `gtg`
        * This generates a set of output files consisting of FASTA / BED / GFF3 
        * With the aaFASTA, we read in each potential sequence and determine if it has a TMD based off of:
            * A min, max range from the start codon (<user_input>).
            * Finding a RegEx of Hydrophobic AAs consecutively next to each other, a total of a <user_input> consecutively together.
                * For instance, if the TMD wanted to be queried for a 12 member hydrophoic region, `-tmdRange = 12`
            * If this fails, it queries to find if there is a Lysine triplet region at the beginning of the desired min range, or at the end (n to n+2 and tmsize-6 to tmsize-2) 
        * OUTPUT : All matching description/sequence pairs are stored, and returned as a `putative_isp.fa` file.
    * `spaninFuncs.py`
        * func `find_lipobox`
        * func `find_tmd`
        * func `tuple_fasta`

## To-dos:
[x] - Determine Methodology for ORF function in `CPT.py`
    * Done. Naive caller used. Table 11 used in functions.

[x] - User input 
    * FASTA (<s>preferred</s>done) ; <s>not entire .gb file. We only want the sequence</s>

[x] - Output putative .fasta files 

[ ] - Mimic returns from `findSpanin.pl` / Best candidate file (?) 
    * These include:
        * Overlap spanins
        * Covered spanins
        * Next-to spanins