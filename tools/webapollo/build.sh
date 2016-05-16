#!/bin/bash
set -ex

APOLLO_HOST=https://cpt.tamu.edu/apollo
APOLLO_USER='hxr@tamu.edu'
APOLLO_PASS=

function annotation_tables(){
    #echo "Fetching Data"
    #python export.py $APOLLO_HOST $APOLLO_USER $APOLLO_PASS $(find student-data/*.gff | sed 's/student-data\///g;s/.gff//g');
    echo "Reporting"
    python ../phage/phage_annotation_table.py \
        --reportTemplateName phageqc_report_annotation_table.tsv \
        --annotationTableCols "rid,name,start,end,strand,length,sd_seq,sd_spacing,start_codon,ig_dist,notes,product,dbxrefs" \
        out.gff3 out.fa > out.tsv;

    echo "Postprocessing"
    awk -F'\t' '{print > "annotation_table_"$1".tsv"}' out.tsv;

    mv annotation_table_\#*.tsv at_header.tsv;

    for i in annotation_table_*.tsv;
    do
        cat at_header.tsv $i > final_$(basename ${i} .tsv)_v2.tsv;
        rm $i;
    done
}

annotation_tables
