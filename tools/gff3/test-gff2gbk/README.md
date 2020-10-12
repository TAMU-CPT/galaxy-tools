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