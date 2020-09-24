# CPT Galaxy Tools

This is the collection of tools developed by the [Center for Phage Technology](https://cpt.tamu.edu/) for their Galaxy instance. The tool set is tailored for the analysis and annotation of bacteriophage genomes, but should also work well for prokaryotic genomes. 

Documentation on the use of common tools and workflows can be found at [https://cpt.tamu.edu/training-material](https://cpt.tamu.edu/training-material).

## Tools

These tools are installed on Galaxy. Some came directly from the [Galaxy Toolshed](). Most are written in, or are compatible with, Python. For general dependencies and tool requirements, see [requirements.txt](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/requirements.txt). Exceptions and special dependencies are described in the wrapper (.xml) for a tool, a few notable cases are listed below.

- The BLAST wrapper uses the [TaxonKit](https://github.com/shenwei356/taxonkit) for TaxID mapping to species level, as supported by BLAST+ versions 2.8.1 and above using the BLASTDBv5 database format
- PAUSE has additional [requirements](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/tools/pause/requirements.txt)
- JBrowse dependencies, including perl, are listed in its [tool_dependencies.xml](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/tools/jbrowse/tool_dependencies.xml)

Tools in the /tools directory are generally organized according to input file type or specific application. Global macros are listed in [macros.xml](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/tools/macros.xml), some additional macros are in relevant sub-directories.

The most commonly used tools have test cases built into the wrappper for code validation via [Planemo](https://github.com/galaxyproject/planemo). Test data for each case is specified, and usually lives in a sub-directory.

Tools in this repository are designed for use inside of a Galaxy instance. Any use or installation independent of Galaxy, or outside the scope detailed for any given tool, is not supported. The CPT Galaxy instance hosts multiple tools produced by external groups which are not present in this repository. Some of these tools require additional permissions or licenses; please contact the relevant tool author(s) if you wish to host these on your Galaxy instance.

#### Deprecation Candidates

Some tools are being set for deprecation and may be removed in future updates.

Tool Name | Script Location | Reason
----------|-----------------|--------
EFetch | edirect/efetch.xml | To be replaced with IUC version

## Testing and Test Data
All tools present in the Structural, Functional, and Comparative workflows are testable with the `planemo test` command or via Galaxy's built in `run_tests.sh` script using the data in each subdirectory's `test-data` folder. The following tools are exceptions to this, as they rely on externally installed tools or services that are outside of Galaxy or Conda installed Packages. 

We are currently working to include these tests to all tools.

Tool Name | Script Location
----------|-----------------
TMHMM (GFF3) | external/tmhmm.xml
[5.33-72.0] Interproscan functional predictions of ORFs | external/interproscan-5.33.xml
JBrowse 0.6.3+cpt | jbrowse/jbrowse.xml
progressiveMauve | comparative/progressivemauve.xml
All webapollo tools | webapollo/*

## Contribution

When contributing new tools or bugfixes to this repo, ensure that you apply [Black code formatting](https://github.com/psf/black) first.


## LICENSE

Except where otherwise noted, all code is [GNU GPLv3](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/LICENSE).
Exceptions include the [webapollo tools](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/tools/webapollo/README.md) and [JBrowse](https://github.tamu.edu/CPT/Galaxy-Tools/tree/master/tools/jbrowse). Genemark also has its own [license](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/tools/genemark/LICENSE).

## Support and Funding

Funding for the CPT Galaxy project comes from the National Science Foundation (awards DBI-1565146, EF-0949351). Additional support has been provided by an Initial University Multidisciplinary Research Initiative supported by Texas A&M University and Texas AgriLife Research, and from the Department of Biochemistry and Biophysics at Texas A&M University. Visit the [CPT website](https://cpt.tamu.edu/people/) for contact information of the Principal Investigators on this project, Ry Young and Jason Gill. 
