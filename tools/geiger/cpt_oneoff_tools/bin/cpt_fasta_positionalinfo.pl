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

use CPT::GalaxyGetOpt;
use Bio::Tools::SeqStats;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file|f' => 'Input file',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['fasta'],
			}
		],
		[
			"amino_acid" => "Amino acid for positional information",
			{
				validate => 'String',
				required => 1,
				default => 'C',
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
				default        => 'posinfo',
				data_format    => 'text/tabular',
				default_format => 'CSV'
			}
		],
	],
	'defaults' => [
		'appid'   => 'AAPosInfo',
		'appname' => 'Amino Acid Positional Information',
		'appdesc' => 'for an input GBK file',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/multi.cds.fa',
			},
			outputs => {
				'results' => ['posinfo.data.csv', 'test-data/outputs/posinfo.default.csv' ],
			}
		},
		{
			test_name    => "Normalized",
			params => {
				'file' => 'test-data/inputs/multi.cds.fa',
				'amino_acid' => 'C',
				'normalize' => '',
			},
			outputs => {
				'results' => ['posinfo.data.csv', 'test-data/outputs/posinfo.normalized.csv' ],
			}
		},
	]
);

# Load the fasta
use CPT::Bio;
my $bio = CPT::Bio->new();
my $fasta_seqio = $bio->getSeqIO($options->{file});
my ($seq_id, $seq);
my @data;
while(my $seqobj = $fasta_seqio->next_seq()){
	my @z = ($seqobj->display_id());
	my $seq = $seqobj->seq();
	my @sseq = split(//,$seq);
	for(my $i=0;$i<scalar(@sseq);$i++){
		if($sseq[$i] eq $options->{amino_acid}){
			if($options->{normalize}){
				push(@z, $i+1 / scalar(@sseq));
			}else{
				push(@z, $i+1);
			}
		}
	}
	push(@data, \@z);
}

my %results = (
	'data' => {
		'header' => [qw(SeqID Data)],
		'data'   => \@data
	},
);

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
        name => 'results',
        GGO => $ggo,
);
$output->CRR(data => \%results);

=head1 DESCRIPTION

This tool produces a table of each position in a sequence a given amino acid can be found . This data can be automatically normalized if that is an appropriate transform for your data.

=cut
