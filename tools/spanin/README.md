# Overview

## The goal of this(these) tool(s) is to take putative protein sequences and determine if any of them are quality candidates for a potential spanin pair (OSP/ISP). 
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
* `generate-putative-osp.py`
    * INPUT : Genomic FASTA
    * METHODOLOGY 
        * Uses the `OrfFinder` object and locates __ALL__ potential start sequences, based on `ttg` / `atg` / `gtg`
        * This generates a set of output files consisting of FASTA / BED / GFF3 
        * With the aaFASTA, we read in each potential sequence and determine if it has a lipobox based off of: 
            * Min distance from start codon to lipobox
            * Max distance from start codon to lipobox
        * Within that range, a regular expression (RegEx) is used to determine if a lipobox is found.
    * OUTPUT : All matching description/sequence pairs are stored, and returned as a `putative_osp.fa` file.
* `generate-putative-isp.py`
    * INPUT : Genomic FASTA
    * METHODOLOGY
    * Uses the `OrfFinder` object and locates __ALL__ potential start sequences, based on `ttg` / `atg` / `gtg`
    * This generates a set of output files consisting of FASTA / BED / GFF3 
    * With the aaFASTA, we read in each potential sequence and determine if it has a TMD based off of:
        * Searches for Lysine snorkels (see `find_tmd`)
        * Searches for large spanning hydrophobic regions (see `find_tmd`)
    * OUTPUT : All matching description/sequence pairs are stored, and returned as a `putative_isp.fa` file.
* `spaninFuncs.py`
    * func `find_lipobox`
        * Uses a choice of (currently two) regular expressions that find lipoboxes upstream of the input sequence. If this passes, the description/sequence is saved.
    * func `find_tmd`
        * Does two primary things:
            1. Ciphs through looking for Lysine snorkels, of a range 1 - 7 upstream of the selected sequence. If it finds a lysine, it looks to see if there is a hydrophobic neighbor, and then looks upstream for a range of transmembrane size - 6 to transmembrane size - 1, and verifies if there is another snorkeling lysine with a hydrophobic neighbor. If this passes, the description/sequence is saved.
            2. Looks for a repeated hydrophobic region within the sequence, based on a range of user inputs (example, transmembrane size 6 - 20, will look for [hydrophobic-AAs]x6 through [hydrophobic-AAs]x20). If this passes, the description/sequence is saved.
    * func `tuple_fasta`
        * Outputs a tuple which contains the description header from the original fasta generation with candidate sequences.

## To-dos:
[x] - Determine Methodology for ORF function in `CPT.py`
    * Done. Naive caller used. Table 11 used in functions. Remove 'I' starts within `tuple_fasta`

[x] - User input 
    * FASTA (<s>preferred</s>done) ; <s>not entire .gb file. We only want the sequence</s>

[x] - Output (correct) putative .fasta files containing ISP/OSP for lambda & T7
    * Works, where for both, the sequence can be found for each of their isp and osp.

[ ] - Mimic returns from `findSpanin.pl` / Best candidate file (?) 
    * These include:
        * Overlap spanins
        * Covered spanins
        * Next-to spanins