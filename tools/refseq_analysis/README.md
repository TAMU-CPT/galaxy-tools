# Summary
Implements a set of piped calls to NCBI using Entrez to confirm if a given input accession has a RefSeq record.

# Process
**Input**: New line file with separated genome accessions _and/or_ gids.

**Methods**: 
* Quality control over the input list, such as white spacing and other basic checks.
* Essentially it follows this command line input:
```bash
esearch -db nuccore -query 'KC821634.1 FR687252.1' | elink -related -name nuccore_nuccore_gbrs | esummary | xtract -pattern
```
* Instead of bash, I wrote this with python using biopython so that it can be utilized in a collapsed package (down the line).

**Output**: Tab delimited table reporting if an input accession has a refseq record and the name of the organism (sanity check). 