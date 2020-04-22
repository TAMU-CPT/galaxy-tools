# CPT Galaxy Tools

This is the collection of tools found in the [Center for Phage Technology](https://cpt.tamu.edu/) Galaxy instance. The entire set is tailored for phage genome annotation and prokaryotic host analysis. 

Documentation for common-use tools can be found at [https://cpt.tamu.edu/training-material](https://cpt.tamu.edu/training-material).

## Tools

These tools are installed on Galaxy. Some came directly from the [Galaxy Toolshed](). Most are written in/compatible with Python. For general required dependencies, see [requirements.txt](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/requirements.txt). Exceptions will be noted in the wrapper (.xml) for a tool, a few cases are listed below.

- BLAST wrapper uses the TaxonKit for TaxID mapping to species level
- PAUSE has additional [requirements](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/tools/pause/requirements.txt)
- JBrowse dependencies, including perl, in its [tool_dependencies.xml](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/tools/jbrowse/tool_dependencies.xml)

Tools in the /tools directory are generally organized according to input file type or specific application. Global macros are listed in [macros.xml](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/tools/macros.xml), some additional macros are in relevant sub-directories.

The most commonly used tools have test cases built into the wrappper for code validation via [Planemo](https://github.com/galaxyproject/planemo). Test data for each case is specified, and usually lives in a sub-directory.

Tools in this repository are meant to be used on Galaxy. Any use or installation independent of Galaxy, or outside the scope detailed for any given tool, is not guaranteed. 


## Contribution

When contributing new tools or bugfixes to this repo, ensure that you apply [Black code formatting](https://github.com/psf/black) first.


## LICENSE

Except where otherwise noted, all code is [GNU GPLv3](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/LICENSE).
Exceptions include the [webapollo tools](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/tools/webapollo/README.md) and [JBrowse](https://github.tamu.edu/CPT/Galaxy-Tools/tree/master/tools/jbrowse). Genemark has its own [license](https://github.tamu.edu/CPT/Galaxy-Tools/blob/master/tools/genemark/LICENSE).

## Support and Funding

Funding for the CPT Galaxy project comes from the National Science Foundation. Additional suuport has come from an Initial University Multidisciplinary Research Initiative supported by Texas A&M University and Texas AgriLife, and from the Department of Biochemistry and Biophysics at Texas A&M University
