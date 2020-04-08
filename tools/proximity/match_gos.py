##### Because I messed up the output, and I don't want JRR to redo her parsing of the output go-synonyms. I am going to store the GO terms
##### that she selected in the originhal "fetch". And then will rerun goQuery and parse out the same lines that she did not want. The new
##### output grom goQuery will retain the keys from the insert database, thus making the mapping FAR less of a hastle.

import pandas as pd

jrr_file = "go-synonym-results-JRR.txt"
og_file = "go-synonym-results.txt"

cols = ["query_term", "GO:id", "GO:name", "description", "restriction", "synonyms"]
jrr_frame = pd.read_csv(jrr_file, names=cols, sep="\t")
print(jrr_frame.head(5))
gos_to_grab = jrr_frame["GO:id"].tolist()
for go in gos_to_grab:
    print(go)