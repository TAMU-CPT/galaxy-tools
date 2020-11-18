PAUSE
=====

PAUSE stands for Pile-up Analysis Using Starts & Ends. It leverages high
throughput read data to determine boundaries of high read coverage
regions with base level granularity. This enables the accurate
determination of the locations of Short Terminal Repeats. Additionally,
we have been able to use this method to accurately position Long
Terminal Repeats, a feat previously much harder to accomplish.

Requirements
------------

-  `Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`__
-  `SAM Tools <http://samtools.sourceforge.net/>`__
-  `R <http://www.r-project.org/>`__
-  `libcpt <https://cpt.tamu.edu/gitlab/cpt/libcpt.git>`__
-  Perl Modules:
-  Bio::DB::Sam Bio::SeqIO CPT CPT::OutputFiles CPT::Report::HTML
   Data::Dumper ExtUtils::MakeMaker File::ShareDir
   File::ShareDir::Install File::Spec File::Temp IPC::Run3
   List::MoreUtils List::Util Moose POSIX Statistics::Descriptive SVG
   Test::More

How to use PAUSE
----------------

Run ``perl pause_analysis.pl --help`` for a full list of options. To run
with the example data:

.. code:: bash

    perl bin/pause_analysis.pl --bam_sense t/Angus.sense.sam.bam.sorted.bam --bai_sense t/Angus.sense.sam.bam.sorted.bam.bai --bam_antisense t/Angus.antisense.sam.bam.sorted.bam --bai_antisense t/Angus.sense.sam.bam.sorted.bam.bai --genome t/Angus.fa

This produces an HTML file (``pause.html``) and a report on the list of
peaks (``pause-extra.csv``)

The PAUSE results consist of two columns, the left column showing a
residue by residue histogram of starts, and the right side showing a
small image of the surrounding area with coverage density and start
peaks labelled.

Once you've done this, you can read `the Analysis Guide <ANALYSIS.md>`__
for further information on how to understand your data.

PAUSE, in depth
---------------

PAUSE simply looks at your mapped data and extracts the starts, ends,
and coverage density from it. With this information in hand, it runs
Continuous Wavelet Transform on the data, in order to identify peaks.
These peaks are then plotted for visual convenience and displayed to the
user.

Some additional programmes are provided to generate standalone plots if
that is needed (``read-sam.pl`` and ``read-sam-zoom.pl``)

PAUSE Suite
-----------

PAUSE consists of several different Perl utilities, found in the
``bin/`` directory:

-  ``bin/pause_analysis.pl``: run the actual PAUSE analysis
-  ``bin/pause_histo.pl``: generates a histogram of the number of reads
   mapping to each base of the genome
-  ``bin/pause_read-sam.pl``: static BAM alignment visualisation. This
   functions much like IGV in that it can show you the # of reads mapped
   to each base in a clean and easy to read way.
-  ``bin/pause_read-sam-zoom.pl``: generates zoomed images of the above

Data Preparation
----------------

Preparing data for PAUSE is relatively trivial. We have galaxy workflows
on hand for anyone using galaxy to ease the process.

.. code:: bash

    # Build the bowtie index
    bowtie2-build genome.fa genome;
    # Run the mapping
    bowtie2 -p 4 -x genome -1 unaligned_1.dat -2 unaligned_2.dat -I 0 -X 250 --very-sensitive --gbar 4 > mapped.bam

    # Filter
    samtools view -S -h -x -f 0 -F 4 -q 0 mapped.bam > mapped.filtered.bam

    # Separate out stuff which maps to each strand
    samtools view -S -b -h -x -f 0 -F 16 -q 0 mapped.filtered.bam > mapped.sense.bam
    samtools view -S -b -h -x -f 16 -F 0 -q 0 mapped.filtered.bam > mapped.antisense.bam

    # Sort it
    samtools sort mapped.sense.bam mapped.sense.sorted.bam
    samtools sort mapped.antisense.bam mapped.antisense.sorted.bam

    # and index it
    samtools index mapped.sense.sorted.bam
    samtools index mapped.antisense.sorted.bam

    # Run PAUSE, which generates a static mapping of the reads to the genome
    perl pause_read-sam.pl --fasta genome.fa --sense_bam mapped.sense.sorted.bam --antisense_bam mapped.antisense.sorted.bam \
        --kb_per_row 20 --line_height 200 --plot_width 1500 --x_border 100 --y_border 150 --plot 8869.dat \
        --plot_files_path . --plot_format HTML

    # Run PAUSE:Analaysis, which generates the table of peaks.
    perl pause_analysis.pl --bam_sense mapped.sense.bam --bam_antisense mapped.antisense.bam \
        --bai_sense mapped.sense.bam.bai --bai_antisense 8868.dat --genome genome.fa \
        --starts_threshold 10

