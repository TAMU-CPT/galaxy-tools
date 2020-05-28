# Overview
This is a modified version of the efetch tool from NCBI and Galaxy. It is a reduced version of efetch, in that it limits the databases and return types. This can be expanded in the future, however, it is currently (05.14.2020) _just_ able to access the protein and nucleotide database. Genbank and Fasta files are the only return file types. 

The power in this tool is in fetching large amounts of files. There is built in sleep functions that will delay large queries to NCBI, as well as attempt to resubmit GET requests if a HTTP error occurs.

Due to not _wanting_ to cause too much of a ruccus if multiple user's are wanting to use this tool, I believe that in our Galaxy config we should limit only one concurrent use of this tool. It is capable of being abused, and since NCBI __will__ block by IP, I think it's worth our time and effort to ensure we do not overwhelm their systems as well as not be on the recieving end of a ban. Helena from the Galaxy team has their efetch tools with the following Galaxy config:

``` xml
<destination id="entrez" runner="local">
</destination>
<limit type="concurrent_jobs" id="entrez">1</limit>
<tools>
  <tool id="ncbi.eutils.efetch" destination="entrez" />
  <tool id="ncbi.eutils.esearch" destination="entrez" />
  <tool id="ncbi.eutils.epost" destination="entrez" />
  <tool id="ncbi.eutils.elink" destination="entrez" />
  <tool id="ncbi.eutils.einfo" destination="entrez" />
  <tool id="ncbi.eutils.esummary" destination="entrez" />
</tools>
```

# Current Implmentation
* Databases: `protein`, `nucleotide`
* Return Types: `fasta`, `genbank`
* Output formats: `fasta`, `genbank`, `multifasta`, `multigenbank`