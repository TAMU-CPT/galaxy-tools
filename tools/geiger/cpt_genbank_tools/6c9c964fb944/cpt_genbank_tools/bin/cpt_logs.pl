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

use CPT;
use Data::Dumper;

my $libCPT  = CPT->new();
my $options = $libCPT->getOptions(
	'options' => [
		[ 'file|f', 'Input file',
			{
				validate => 'File/Input',
				#file_format => ['Genbank'],
				required => 1,
			}
		],
	],
	'outputs' => [
		[
			'results',
			'LOGS Results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'logs',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'LOGS',
		'appname' => 'LOGS',
		'appdesc' => 'calculates length of genomic sequences. Useful for normalising genomic data.',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				"file" => "test-data/inputs/multi.gbk",
				"results_format" => "CSV",
			},
			outputs => {
				"results" => ["logs.Sheet1.csv", "test-data/outputs/logs.csv"]
			},
		},
	],
);

# Which tags should we consider as part of PIGS? These should be coding sequences
my %wanted_tags = (
	CDS => 1,
);

my %data = (
	'Sheet1' => {
		header => ['Sequence', 'Length'],
		data => [],
	}
);
my @data_data;

use Bio::SeqIO;
use CPT::Bio;
my $bio = CPT::Bio->new();
my $seqio_object = $bio->getSeqIO($options->{file});
# For all genomes in the GBK file
while(my $seq_object = $seqio_object->next_seq){
	push(@data_data, [ $seq_object->display_id(), $seq_object->length() ] );
}

$data{Sheet1}{data} = \@data_data;
use CPT::OutputFiles;
my $csv_output = CPT::OutputFiles->new(
        name => 'results',
        libCPT => $libCPT,
);
$csv_output->CRR(data => \%data);

=head1 NAME

Length of Genomic Sequences (LOGS)

=head1 DESCRIPTION

Exactly what it says on the tin. Calculates length of all genomes in a genbank file. This is mostly useful for pipelines

=cut
