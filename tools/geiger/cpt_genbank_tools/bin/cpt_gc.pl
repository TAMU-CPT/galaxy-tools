#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[ 'file|f', 'Input file',
			{
				validate => 'File/Input',
				file_format => ['genbank', 'embl', 'txt'],
				required => 1,
			}
		],
		[ 'tag', 'Genomic tag of interest to analyse',
			{
				validate => 'Genomic/Tag',
				multiple => 1,
			}
		],
	],
	'outputs' => [
		[
			'results',
			'GC% Results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'gc',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'Percent_GC',
		'appname' => 'GC Percentage',
		'appdesc' => 'calculates %%GC in a genome, as well as across all sequences in the genome',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default GBK",
			params => {
				'file' => 'test-data/inputs/multi.gbk',
			},
			outputs => {
				'results' => [ 'gc.Sheet1.csv', 'test-data/outputs/gc.default.csv'],
			}
		},
		{
			test_name    => "Default GBK CDSs",
			params => {
				'file' => 'test-data/inputs/multi.gbk',
				'tag' => 'CDS',
			},
			outputs => {
				'results' => [ 'gc.Sheet1.csv', 'test-data/outputs/gc.cds.csv'],
			}
		},
		{
			test_name    => "Default .fa",
			params => {
				'file' => 'test-data/inputs/multi.cds.fa',
			},
			outputs => {
				'results' => [ 'gc.Sheet1.csv', 'test-data/outputs/gc.fa.csv'],
			}
		},
	],
);

# Which tags should we consider as part of PIGS? These should be coding sequences
my %data = (
	'Sheet1' => {
		header => ['Parent Sequence', 'Sequence ID', 'GC %'],
		data => [],
	}
);
my @data_data;

my %wanted_tags = map {$_ => 1} @{$options->{tag}};

use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;
my $seqio_object = $bio->getSeqIO($options->{file});
# For all genomes in the GBK file
while(my $seq_object = $seqio_object->next_seq){
	push(@data_data, [ $seq_object->display_id(), 'Whole Genome', analyse_sequence($seq_object->seq)]);
	for my $feat_object ($seq_object->get_SeqFeatures) {
		if($wanted_tags{ $feat_object->primary_tag }){
			push(@data_data, [ $seq_object->display_id(), $bio->_getIdentifier($feat_object), analyse_feature($feat_object)]);
		}
	}
}

$data{Sheet1}{data} = \@data_data;
use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
        name => 'results',
        GGO => $ggo,
);
$output->CRR(data => \%data);

sub analyse_feature {
	my ($feat_object) = @_;
	my $loc = $feat_object->location;
	if(ref($loc) eq 'Bio::Location::Simple'){
		return analyse_sequence($feat_object->seq->seq());
	}else{
		return analyse_sequence($feat_object->spliced_seq->seq());
	}
}
sub analyse_sequence {
	my ($sequence) = @_;
	my %bases;
	foreach my $char(split //,$sequence){
		$bases{$char}++;
	}
	return ($bases{C} + $bases{G}) / ($bases{A} + $bases{C} + $bases{G} + $bases{T} );
}

=head1 NAME

% GC Calculator

=head1 DESCRIPTION

This calculator will determine % GC across the whole genome, as well as for every feature in the genome.

You can read more about calculating % GC on L<wikipedia|https://en.wikipedia.org/wiki/GC-content>

=cut
