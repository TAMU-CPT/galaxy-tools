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

my $libCPT = CPT->new();

my $options = $libCPT->getOptions(
	'options' => [
		[
			'file' => 'Input file',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
		[
			"tag" => "Analyse by a
			specific tag",
			{
				validate => 'Genomic/Tag',
				required => 1,
				default  => 'CDS',
			}
		],
		[ 'normalize', 'Normalize all indices over the length of the sequence', { validate => 'Flag' } ],
	],
	'outputs' => [
		[
			'report',
			'First index of AAs',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'report',
				data_format    => 'text/tabular',
				default_format => 'CSV'
			}
		],
	],
	'defaults' => [
		'appid'   => 'FirstIndex',
		'appname' => 'First Index of AAs',
		'appvers' => '1.94',
		'appdesc' => 'for an input GBK file',
	],
	'tests' => [
		{
			test_name        => "Default from GBK",
			params           => {
				'file'   => 'test-data/inputs/multi.gbk',
				'tag'    => 'CDS',
			},
			outputs          => {
				'report' => ['report.Index.csv', 'test-data/outputs/firstindex.default.csv'],
			}
		},
		{
			test_name        => "Default from fasta",
			params           => {
				'file'   => 'test-data/inputs/multi.cds.fa',
				'tag'    => 'CDS',
			},
			outputs          => {
				'report' => ['report.Index.csv', 'test-data/outputs/firstindex.fa.csv'],
			}
		},
		{
			test_name        => 'Normalized from gbk',
			params           => {
				'file'   => 'test-data/inputs/multi.gbk',
				'tag'    => 'CDS',
				'normalize'=>'',
			},
			outputs          => {
				'report' => ['report.Index.csv', 'test-data/outputs/firstindex.norm.csv'],
			}
		},
	],
);

my %args = (
	'file'      => $options->file,
	'callback'  => \&func,
	'translate' => 1,
	'header'    => 1,
	'subset'    => $options->{tag},
);
use CPT::Bio;
my $cptbio = CPT::Bio->new();
$cptbio->parseFile(%args);


sub func {
	my $response_ref = shift;
	my @response     = @{$response_ref};

	use CPT::BioData;
	my $bd    = CPT::BioData->new();
	my %table = %{$bd->getTranslationTable()};
	my %reverse_table;
	foreach ( keys %table ) {
		push( @{ $reverse_table{ $table{$_} } }, $_ );
	}
	my @alphabet = sort(keys(%reverse_table));

	my @header = qw(ID);
	if($options->{normalize}){
		push(@header, 'Length');
	}
	push(@header,@alphabet);

	# prepare our results table. 
	my %results = (
		'Index' => {
			'header' => \@header,
			'data' => [],
		},
	);

	my @data;
	foreach (@response) {
		my ( $header, $sequence ) = @{$_};
		my @line = ($header);

		if($options->{normalize}){
			push(@line, length($sequence));
		}

		foreach my $aa (@alphabet) {
			my $idx = index( $sequence, $aa );
			if($options->{normalize}){
				$idx /= length($sequence);
			}
			push( @line, $idx );
		}
		push( @data, \@line );

	}
	$results{'Index'}{'data'} = \@data;

	use CPT::OutputFiles;
	my $output = CPT::OutputFiles->new(
		name => 'report',
		libCPT => $libCPT,
	);
	$output->CRR(data => \%results);
}

=head1 DESCRIPTION

This tool produces a table of the first time it encounters each amino acid (known to the tool) in each sequence. This data can be automatically normalized if that is an appropriate transform for your data.

=cut
