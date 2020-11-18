#!/usr/bin/perl
#
#       Code written by Eric Rasche
#               mailto:rasche.eric@yandex.ru
#               tel:   404.692.2048
#               http://eric.rasche.co.uk
#       for
#               Center for Phage Technology
#

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use File::Temp qw/ tempdir tempfile /;
use Data::Dumper;
use IPC::Run3;
	use Bio::SearchIO;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file',
			'Input file',
			{
				required => 1,
				validate => 'File/Input',
				multiple => 1,
			}
		],
		[
			'window_size',
			'Length of segments that a genome is chopped into before the blasting process',
			{
				validate => 'Int',
				min => 500,
				required => 1,
				default => 1000,
			}
		],
		[
			'step_size',
			'Step size between each new segment of the genome',
			{
				validate => 'Int',
				min => 200,
				required => 1,
				default => 200,
			}
		],
	],
	'outputs' => [
		[
			'results',
			'ANI Table',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'ani',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appname' => 'ANI',
		'appid'   => 'ANI',
		'appdesc' =>
'calculates average nucleotide identity for a set of genomes',
	],
	'tests' => [
		{
			test_name => "Default",
			params => {
				'file' => 'test-data/inputs/ani.fa',
				'results_format' => 'YAML',
			},
			outputs      => {
				results => ['ani.yml', 'test-data/outputs/ani.yml'],
			}
		},
	]
);


# Process GBK (or fasta) to fasta genomes
my @files = @{ $options->file };
#my $cleanup = !defined($options->{noclean});
my $cleanup = 1;

# Create a temp dir for this run, keeps our data cleaner.
my $dir = tempdir('cpt.ani.XXXXXXX',CLEANUP => $cleanup );

my %fd;

my ( $merged_fasta_fh, $merged_fasta_path ) = tempfile(
	'XXXXXXXXX',
	UNLINK => $cleanup,
	DIR    => $dir,
	SUFFIX => '.fa'
);

# Loop over input files
foreach my $file(@files){
	# Parse out data, possibly will return multiple "files"
	use CPT::Bio;
	my $bio = CPT::Bio->new();
	my %args = (
		'file'      => $file,
		'header'    => 1,
		'subset'    => 'whole',
	);

	foreach my $arrayref($bio->parseFile(%args)){
		# For each sequence
		foreach my $subseqref(@{$arrayref}){
			# For each (possible) subsequence (>=1)
			my %local_file_info;

			# Create temp file with fasta sequence)
			my ($header,$seq) = @{$subseqref};
			my ( $fh, $path ) = tempfile(
				'XXXXXXXXX',
				UNLINK => $cleanup,
				DIR    => $dir,
				SUFFIX => '.fa'
			);
			my $fasta_id = substr($header,1);

			# spliting genome in chunks of $window_size and store in a file specific to that genome
			for(my $i = 1;$i < length($seq); $i += $options->{step_size}){
				print $fh sprintf(">%s.%s\n%s\n", $fasta_id, $i, substr($seq,$i,$options->{window_size}));
			}

			# Also store the complete sequence in a "merged" .fa file
			print $merged_fasta_fh "$header\n$seq\n";

			$fd{$fasta_id} = { path => $path, length => length($seq) };
		}
	}
}

# create blast db
my $db_location = create_blast_db($merged_fasta_path);

my %ani_table;
# run blast comparisons, a given genome in FD against ALL OTHERS (in the merged fasta file)
foreach my $query(keys %fd){
	# All blast hits
	my $outfile = blastn($db_location, $fd{$query}{path});
	# Data struct to hold good ones
	
	my %scoring_table;
	my %reductions;
	my %target_reduction;

	my $in = new Bio::SearchIO(-format => 'blastxml',
		   -file   => $outfile);
	while( my $result = $in->next_result ) {
		## $result is a Bio::Search::Result::ResultI compliant object
		while( my $hit = $result->next_hit ) {
			## $hit is a Bio::Search::Hit::HitI compliant object
			while( my $hsp = $hit->next_hsp ) {
				## $hsp is a Bio::Search::HSP::HSPI compliant object
				for(my $i=$hsp->start('hit'); $i <= $hsp->end('hit'); $i++){
					$target_reduction{$query}{$hit->name}{$i} = 1;
				}
				if( $hsp->length('total') > ($options->{window_size} * .7) && $hsp->percent_identity >= 75 ) {
					# For short sequences, correcting here will drop self-self matches to 85%, uncorrected it's 99%
					$scoring_table{$query}{$hit->name} += $hsp->percent_identity;
					$reductions{$query}{$hit->name}++;
				}
			}
		}
	}
	
	foreach my $b(keys($scoring_table{$query})){
		my $hit;
		foreach my $var(keys($target_reduction{$query}{$b})){
			$hit++;
		}
		my $length = $fd{$query}{length} + $fd{$b}{length};

		my $uncorrected_ani = $scoring_table{$query}{$b} / $reductions{$query}{$b};
		my $correction_factor = 2*$hit/$length;
		$ani_table{$query}{$b} = $uncorrected_ani * $correction_factor;
	}
}


# Reformat results
my @table_data;
my @keys = sort(keys(%fd));
foreach my $query(@keys){
	my @row = ($query);
	foreach my $subject(@keys){
		push(@row, $ani_table{$query}{$subject});
	}
	push(@table_data, \@row);
}
my %response = (
	'PercentageConservedDNA' => {
		header => ['', @keys],
		data => \@table_data,
	}
);

# Output them
use CPT::OutputFiles;
my $mist_output = CPT::OutputFiles->new(
	name => 'results',
	GGO => $ggo,
);
$mist_output->CRR(data => \%response);

sub blastn {
	my ($blastdb_path, $file) = @_;
	my ($in, $out, $err);

	my ( $blastxml_out_fh, $blastxml_out_path ) = tempfile(
		'XXXXXXXXX',
		UNLINK => $cleanup,
		DIR    => $dir,
		SUFFIX => '.xml'
	);

	my @command = ('blastall',
		'-p', 'blastn',
		'-d', $blastdb_path,,
		'-i', $file,
		'-F', 'F',
		'-e', '0.001',
		'-v', '1',
		#'-b', '1',
		# This flag apparently limits number of hits we'll get back
		'-X', '150',
		'-q', '-1',
		'-a', '2',
		'-m', '7',
		'-o', $blastxml_out_path,
	);
	if($options->{verbose}){
		print STDERR join ' ', @command,"\n";
	}
	#my @command = ('blastn',
		#'-xdrop_gap_final', '150',
		#'-gapopen', '2',
		#'-gapextend', '1',
		#'-penalty', '-1',
		#'-db', $blastdb_path,
		#'-query', $file,
		#'-outfmt','6',
		#'-perc_identity', '30'
	#);
	run3 \@command, \$in, \$out, \$err;
	if($err){
		print STDERR $err;
	}
	return $blastxml_out_path;
}

sub create_blast_db{
	my ($file) = @_;
	my ($in, $out, $err);
	my ( $blastdb_fh, $blastdb_path ) = tempfile(
		'XXXXXXXXX',
		UNLINK => $cleanup,
		DIR    => $dir,
		SUFFIX => ''
	);

	my @command = ('makeblastdb', '-dbtype', 'nucl', '-out', $blastdb_path, '-in', $file);
	#my @command = ('formatdb', '-i', $file, '-n', $blastdb_path);
	run3 \@command, \$in, \$out, \$err;
	if($err){
		print STDERR $err;
	}
	return $blastdb_path;
}

=head1 Why You Should Use This and NOTHING Else

Methods like EMBOSS Stretcher (and other sequence alignment) are good for aligning sequences, they are B<NOT> good for comparing the similarity of genomes, as they are B<very> susceptible to rearrangements of the genome, producing many false negatives.

Methods like the previous ANI script fail in that they will report falsely high scores for a small portion of the genomes matching. If you have an aligned region of 1kb in the old version, and 50kb in the new version that are identical, but the other 99% of the genome is I<completely 100% different>, they will report the genomes as being identical.

=head1 NAME

ANI - Average Nucleotide Identity

=head1 PURPOSE

This tool provides a means by which you can compare N-by-N genomes for their overall average nucleotide similarity, avoiding the pitfalls of using Emboss Stretcher/MegaBlast or similar methods.

This tool does an N-by-N comparison of genomes for average nucleotide identity. This tool was based upon previous work of <Konstantinidis and Tiedje|https://www.ncbi.nlm.nih.gov/pmc/articles/PMC549018/>, wherein they described a method by which two genomes could be analysed for their similarity. This method split the genome into subsets and blasted them, before comparing and averaging the score of the "good fragments".

We've modified their algorithm as we believe it did not accurately reflect the biology, and their algorithm was flawed in how they determined the final ANI#. Our tool more accurately reflects the true % similarity.

=head1 How it Works

Each genome is split into subsets, defined by the window and step size parameters. These fragments are blasted against one another, and then scored. Similar fragments (Percent identity > 75 and length of hit > 70% of gene), are collected. These scores are averaged, and corrected for the total % of the target genome they cover.


=cut
