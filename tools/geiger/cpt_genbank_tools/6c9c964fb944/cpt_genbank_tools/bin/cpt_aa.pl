#!/usr/bin/perl
#
#       Code written by Eric Rasche
#               mailto:rasche.eric@yandex.ru
#               tel:404.692.2048
#               http://eric.rasche.co.uk
#       for
#               Center for Phage Technology
#

use strict;
use warnings;

use CPT;
use Bio::Tools::SeqStats;
use CPT::BioData;
use CPT::Bio;

my $libCPT = CPT->new();

my $options = $libCPT->getOptions(
	'options' => [
		[
			'file' => 'Input Genbank file',
			{
				required => 1,
				validate => 'File/Input',
				multiple => 1,
			}
		],
		[ 'normalize' => 'Normalize all data' ],
	],
	'outputs' => [
		[
			'results',
			'Analysis results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'aa',
				data_format    => 'text/tabular',
				default_format => 'CSV'
			}
		],
	],
	'defaults' => [
		'appid'   => 'AA',
		'appname' => 'Amino Acid Statistics',
		'appdesc' => 'for an input GBK file',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/single.gbk',
				'results_format' => 'YAML',
			},
			outputs      => {
				'results' => ['aa.yml', 'test-data/outputs/aa.yml'],
			}
		},
		{
			test_name => "normalized+individual",
			params => {
				'file' => 'test-data/inputs/single.gbk',
				'results_format' => 'YAML',
				'normalize' => " ",
			},
			outputs => {
				'results' => ['aa.yml', 'test-data/outputs/aa-norm-indiv.yml'],
			}
		},
	],
);


# Preparing table header
my @header = ('Sequence');
my @aa_header = ('Sequence');
if($options->{normalize}){
	push(@header, 'Count');
	push(@aa_header, 'Count');
}

# Used in table header
my $bd    = CPT::BioData->new();
my $bio   = CPT::Bio->new();
my %table = %{$bd->getTranslationTable()};

# Keys for header and lookup in SeqStats
my @bases = ();
foreach ( sort(keys(%table)) ) {
	# Store for use in an identically ordered array
	push(@bases, $_);
	# Also store as a formatted header
	push(@header, sprintf("%s (%s)", $_, $table{$_}));
	
} # gets XXX => X

# Keys for AA headers
my @aakeys = ();
my %rt;
foreach ( sort(keys(%table)) ) {
	$rt{$table{$_}} = 1;
}
@aakeys = sort(keys(%rt));
my %short_long_lookup = (
	'G',  => 'Gly',
	'P', => 'Pro',
	'A', => 'Ala',
	'V', => 'Val',
	'L', => 'Leu',
	'I', => 'Ile',
	'M', => 'Met',
	'C', => 'Cys',
	'F', => 'Phe',
	'Y', => 'Tyr',
	'W', => 'Trp',
	'H', => 'His',
	'K', => 'Lys',
	'R', => 'Arg',
	'Q', => 'Gln',
	'N', => 'Asn',
	'E', => 'Glu',
	'D', => 'Asp',
	'S', => 'Ser',
	'T', => 'Thr',
	'X', => 'XXX',
	'*' => 'Stop',
);
push(@aa_header, map { sprintf("%s (%s)", $short_long_lookup{$_}, $_) } @aakeys);

# Final output table
my %data = (
	'CodonUsageInCDSs' => {
		header => \@header,
		#data => [],
	},
	'AAUsageinCDS' => {
		header => \@aa_header,
		#data => [],
	},
	'OverallAAUsage' => {
		header => \@aa_header,
		#data => [],
	},
);

use Bio::SeqIO;
my @codon_usage_data;
my @aa_usage_data;
my %overall_data;
foreach my $file(@{$options->{file}}){
	# For all genomes in the GBK file
	my $results_ref = $bio->parseFile(
		'file'      => $file,
		'header'    => 1,
		'subset'    => ["CDS"],
		'translate' => 0,
	);

	# Foreach CDS in the genome
	foreach my $item(@{$results_ref}){
		my ($header, $sequence) = @{$item};
		# Convert to PrimarySeq for easy use in SeqStats.
		my $seqobj = Bio::PrimarySeq->new(
			-seq      => $sequence,
			-alphabet => 'dna',
			-id       => 'test'
		);
		# Codon counts
		my $codons = Bio::Tools::SeqStats->count_codons($seqobj);
		my @row = ($header);
		my @aa_row = ($header);
		my %aa_row_t;
		# Divide length by three due to NTs vs AAs
		my $length = length($sequence) / 3;
		if($options->{normalize}){
			push(@row, $length);
		}
		foreach my $base (@bases){
			if(defined $codons->{$base}){
				$aa_row_t{$table{$base}} += $codons->{$base};
				if($options->{normalize}){
					push(@row, $codons->{$base} / $length);
				}else{
					push(@row, $codons->{$base});
				}
			}else{
				push(@row, 0);
			}
		}
		# Add row to cud
		push(@codon_usage_data, \@row);
		# Transform @aa_row_t to @aa_row
		foreach(@aakeys){
			if(defined $aa_row_t{$_}){
				push(@aa_row, $aa_row_t{$_});
			}else{
				push(@aa_row, 0);
			}
		}
		push(@aa_usage_data,\@aa_row);

		# Transform @aa_row_t to %overall_data
		foreach(@aakeys){
			if(defined $aa_row_t{$_}){
				$overall_data{$_} += $aa_row_t{$_};
			}
		}
	}
}
$data{'CodonUsageInCDSs'}{'data'} = \@codon_usage_data;
$data{'AAUsageinCDS'}{'data'} = \@aa_usage_data;

use Data::Dumper;
#print Dumper \@aa_usage_data;

my @overall_data_t = ('Sum');
foreach(@aakeys){
	if(defined $overall_data{$_}){
		push(@overall_data_t, $overall_data{$_});
	}else{
		push(@overall_data_t, 0);
	}
}
$data{'OverallAAUsage'}{'data'} = [\@overall_data_t];

use CPT::OutputFiles;
my $csv_output = CPT::OutputFiles->new(
        name => 'results',
        libCPT => $libCPT,
);
$csv_output->CRR(data => \%data);


=head1 NAME

Amino Acid Statistics

=head1 DESCRIPTION

This tool calculate the following statistics:

=over 4

=item Codon usage in CDSs

How often each codon (e.g., AAT) was used in each CDS

=item AA Usage in CDSs

How often each AA residue was used in a given CDS

=item Overall AA Usage

Amino acids used in every CDS (summation of previous table).

=back

=cut
