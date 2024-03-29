# Overview

## The goal of these tools is to take protein sequences from naively called ORFs and determine if any of them are quality candidates for a potential spanin gene (USP) or gene pair (OSP/ISP). 
### Current Scripts:
* `generate-putative-isp.py`
    * INPUT : Genomic FASTA 
    * OUTPUT : Putative candidates for isp
* `generate-putative-osp.py`
    * INPUT : Genomic FASTA
    * OUTPUT : Putative candidates for osp
* `generate-putative-usp.py`
    * INPUT : Genomic FASTA
    * OUTPUT : Putative candidates for usp
* `spaninFuncs.py`
    * Functions that drive the scripts within this directory

## Requirements:
* regex (as of 2.17.2020 not needed...)
* cpt.py (copied and edited for handling specific task for these tools)

## Script Descriptions / Methodologies:
_aside: Much of this code base from isp/osp/usp/findspanin tools is **NOT** very DRY. I would like to go in and clean this up; however it may be improved down the line when we incorporate ML into the pipeline_
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
* `generate-putative-usp.py`
    * INPUT : Genomic FASTA
    * Does a combination of the spanin func lipobox and tmd functions with different default parameters than their isp/osp counterparts.
* `spaninFuncs.py`
    * func `find_lipobox`
        * Uses a choice of (currently two) regular expressions that find lipoboxes upstream of the input sequence. If this passes, the description/sequence is saved.
    * func `find_tmd`
        * Does two primary things:
            1. Ciphs through looking for Lysine snorkels, of a range 1 - 7 upstream of the selected sequence. If it finds a lysine, it looks to see if there is a hydrophobic neighbor, and then looks upstream for a range of transmembrane size - 6 to transmembrane size - 1, and verifies if there is another snorkeling lysine with a hydrophobic neighbor. If this passes, the description/sequence is saved.
            2. Looks for a repeated hydrophobic region within the sequence, based on a range of user inputs (example, transmembrane size 6 - 20, will look for [hydrophobic-AAs]x6 through [hydrophobic-AAs]x20). If this passes, the description/sequence is saved.
    * func `tuple_fasta`
        * Outputs a tuple which contains the description header from the original fasta generation with candidate sequences.
    * There are many more functions and I wont list them all here. I believe that the majority of the functions have appropriate descriptions.

## Release Notes
Functional Prediction>Lysis Gene Prediction: ISP candidates, OSP candidates, USP candidates, Find Spanin
These tools are used in conjunction to help find potential spanins. ISP, OSP, and USP candidate tools take an input genomic fasta, where all putative proteins are gathered, based on having a “ATG”, “TTG”, or “GTG” start site. Each potential gene is fed through a set of criteria, including TMD distance from the start site for i-spanins, and lipobox presence for o-spanins, or both for u-spanins. The putative_isp and putative_osp fasta outputs are fed to the Find Spanin tool. The output is broken up based on their class:embedded, overlapping, or separate.
