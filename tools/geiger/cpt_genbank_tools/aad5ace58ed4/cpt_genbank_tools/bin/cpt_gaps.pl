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
				#file_format => ['Genbank'],
				file_format => ['genbank', 'embl', 'txt'],
				required => 1,
			}
		],
		['mode', 'Operational mode',
			{
				validate => 'Option',
				options => {
					'cds_only' => 'Only look at CDS features when determining gap length',
					'any' => 'All features will be used when determining gap locations',
				},
				required => 1,
				default => 'cds_only',
			}
		],
	],
	'outputs' => [
		[
			'results',
			'Gap Results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'gaps',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'GAPS',
		'appname' => 'GAPS',
		'appdesc' => 'lists gaps in genomic sequences',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/multi.gbk',
			},
			outputs => {
				'my_output_data' => ["gaps.Sheet1.csv", 'test-data/outputs/gaps.Sheet1.csv' ],
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
		header => ['Sequence', 'Start', 'End', 'Length'],
		data => [],
	}
);
my @data_data;

use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;
my $seqio_object = $bio->getSeqIO($options->{file});
# For all genomes in the GBK file
while(my $seq_object = $seqio_object->next_seq){
	# Coverage history, we'll count 0s vs 1s after
	my %coverage;
	foreach(1..$seq_object->length()){
		$coverage{$_} = 0;
	}

	# For all of our features, if we want them, bump coverage #s
	for my $feat_object ($seq_object->get_SeqFeatures) {
		if($feat_object->primary_tag ne 'source'
			&&
			(
			($options->{mode} eq 'cds_only' && $feat_object->primary_tag eq 'CDS')
			||
			($options->{mode} eq 'any')
			)
		){
			my $loc = $feat_object->location;
			if(ref($loc) eq 'Bio::Location::Simple'){
				foreach(my $i = $feat_object->start; $i < $feat_object->end; $i++){
					$coverage{$i} = 1;
				}
			}else{
				for my $location ( $loc->sub_Location ) {
					foreach(my $i = $feat_object->start; $i < $feat_object->end; $i++){
						$coverage{$i} = 1;
					}
				}
			}
		}
	}


	my $last_state = $coverage{1};
	my $string='';
	for(my $i = 1; $i<$seq_object->length();$i++){
		$string .= $coverage{$i};
	}

	while($string =~ /([0]+)/g){
		push(@data_data, [$seq_object->display_id(), $-[0], $+[0], abs($-[0] - $+[0])]);
	}
}

$data{Sheet1}{data} = \@data_data;
use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'results',
	GGO => $ggo,
);
$crr_output->CRR(data => \%data);

=head1 NAME

Genomic OverlAp caLculator (GOAL)

=head1 DESCRIPTION

Locate regions of overlap between genes.

=cut
