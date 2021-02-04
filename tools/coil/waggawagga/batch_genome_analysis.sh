#!/usr/bin/env bash
# Make batch SAH-analysis of all genome datasets in waggawagga-cli-0.5.x/test
# cd /fab8/dosi/server/waggawagga-cli-0.5.0-linux-x86_64
# Only FASTA-files with .fa-extension will be processed
for file in test/*.fa;
do 
  # Preparation of filename
  filename=$(basename "$file")
  extension="${filename##*.}"
  filename="${filename%.*}"
  echo "Make SAH-analysis for FASTA-file '$file' "

#  # Complete anaysis with clean-up
#  echo "./waggawagga-cli -a complete -g -f $file --tidy"
#  ./waggawagga-cli -a complete -g -f "${file}" --tidy

#  # Re-analysis with clean-up
#  echo "./waggawagga-cli -a evaluation -i $filename --tidy"
#  ./waggawagga-cli -a evaluation -i "${filename}" --tidy

#  # Re-analysis without clean-up
#  echo "./waggawagga-cli -a evaluation -i $filename"
#  ./waggawagga-cli -a evaluation -i "${filename}"

	# Post analysis of cutoff-filter CSV
	echo "results/select_unique_sahs.rb $file results/$filename/sah_score_14avg_window_14_cutoff_filter.csv"
	results/select_unique_sahs.rb "${file}" results/"${filename}"/sah_score_14avg_window_14_cutoff_filter.csv > results/"${filename}"/post_sah_analysis.txt
done
