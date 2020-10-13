# Explanation
* In handeling the update to the newer Galaxy and Apollo, we needed a way to vet a few new changes:
    * The old RBS "kludge" still works on **old (pre update)** annotated genomes --> `shady_old.fa + shady_old.gff3 = shady_old.gbk`
        * `shady_old_confirm.gbk` is *from* the older run 
    * "Shine_Dalgarno_sequence" features will be caught and converted to RBS --> `sun.fa + sun.gff3 = sun.gbk`
    * Mu, lambda, and T1 each have unique features that are extracted.
        * Features that have been added:
            * "sequence_feature" --> "misc_feature"
            * "recombination_feature" --> "misc_recomb"
            * "sequence_alteration" --> "variation"
            * "binding_site" --> "misc_binding"
* `generate_gbks.sh` will regenerate the genbanks for testing (`bash generate_gbks.sh`)
* Lastly, the conditional in the `handle_non_gene_features` function will extract all extra features to the output genbank. 
    * Mu has an extra feature (see `handle_non_gene_features` function) that is caught by the conditional modification. Test implemented to verify it is added to genbank.