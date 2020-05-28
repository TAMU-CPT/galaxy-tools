# Overview
> input: multi protein fasta file
> output: multi fasta of positive candidates and a table summarizing the stats for each candidate with identifier, length of potential SAR, topology orientation, calculate %G and %A

# Requirements
* python 3.6+
* biopython
* pandas
* numpy

# Outline
1. Read in input multi fasta
    * multi fasta parsed by `biopython_parsing.py`
2. Check SAR requirements
    * Min length check
    * Max length check (user dictated)
    * Hydrophobic residues (Ile, Leu, Val, Phe, Tyr, Trp, Met) except often rich in Gly, Ala, and/or Ser residues
        * <s>Option 1: FIWLVMYAGS</s> Using option #2
        * Option 2: FIWLVMYCATGS # add C and T
    * Lysines can be present in the hydrophobic stretch if within 3 residues of the domain boundary (lysine snorkeling)
    * Topology check
        * N term (net positive charge)
        * C term catalytic domain
3. Return candidates and multi fasta.
4. Write statistics to output file in table format
    * identifier :: length of peptide :: topology orientation :: %G and %A :: likely more later

# Testing
* Mu (mu-proteins.fa) for a TP