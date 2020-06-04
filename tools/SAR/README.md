# Overview
> input: multi protein fasta file
> output: multi fasta of positive candidates and a table summarizing the stats for each candidate with identifier, length of potential SAR, topology orientation, calculate %G and %A

# Requirements
* python 3.6+
* biopython
* <s>pandas</s> --> thought I might use this, but ended up seeing I wouldnt needed it by the time I was wrapping this up.
* <s>numpy</s>

# Outline
1. Read in input multi fasta
    * multi fasta parsed by `biopython_parsing.py`
2. Check SAR requirements
    * <s>Min length check</s> --> Currently omitted
    * <s>Max length check (user dictated)</s> --> Currently omitted
    * Hydrophobic residues (Ile, Leu, Val, Phe, Tyr, Trp, Met) except often rich in Gly, Ala, and/or Ser residues
        * <s>Option 1: FIWLVMYAGS</s> Using option #2
        * __Option 2: FIWLVMYCATGS # add C and T --> This is what is being used__
    * Lysines can be present in the hydrophobic stretch if within 3 residues of the domain boundary (lysine snorkeling)
        * <s>Currently, I'm not checking for lysine snorkels as the hydrophobic region present will still be caught.</s>
        * Snorkelers are found by if a Lys is on the first or last index of the sequence range being inspected, checks for hydrophobic residues between it and either the beginning or end of the sequence.
    * Topology check
        * N term (net positive charge)
        * C term catalytic domain
3. Return candidates and multi fasta.
4. Return candidates in multi gff3.
5. Write statistics to output file in table format
    * identifier :: length of peptide :: topology orientation :: %G and %A :: likely more later
    * Been reworked to include what is currently in tab-separated format:
        * ["Name","Protein Sequence","Protein Length","SAR Length","Putative SAR Sequence","SAR Start Location","[res%]","N-term Sequence","N-term net Charge"]

# File Summaries
* `SAR_functions.py`
    * Has the SAR class and accompanying methods
* `SAR_finder.py`
    * The executed script by Galaxy
* `biopython_parsing.py`
    * _might scale to a_ sym link for parsing bio related files, otherwise will just be related to addressing this experiment.
* `file_operations.py`
    * _mich scale to a_ sym link for operating on files and exporting them, otherwise will just be related to addressing this experiment.

# Testing
* Mu (mu-proteins.fa) for a TP