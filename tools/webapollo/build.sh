APOLLO=https://cpt.atmu.edu hxr@tamu.edu PASSWORD
function annotation_tables(){
    python export.py $APOLLO $(find student-data/*.gff | sed 's/student-data\///g;s/.gff//g');
    python ../phage/phage_annotation_table.py \
        --reportTemplateName phageqc_report_annotation_table.tsv \
        --annotationTableCols "rid,name,start,end,strand,length,sd_seq,sd_spacing,start_codon,ig_dist,notes,product,dbxrefs" \
        out.gff3 out.fa > out.tsv;

    awk -F'\t' '{print > "annotation_table_"$1".tsv"}' out.tsv;

    mv annotation_table_\#*.tsv at_header.tsv;

    for i in annotation_table_*.tsv;
    do
        cat at_header.tsv $i > final_$i;
    done
}
