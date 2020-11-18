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
	],
	'outputs' => [
		[
			'results',
			'GOAL Results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'goal',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'GOAL',
		'appname' => 'GOAL',
		'appdesc' => 'calculates overlap of annotated sequences in a genome',
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
				"results" => ["goal.Sheet1.csv", "test-data/outputs/goal.csv"]
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
		header => ['Hit from', 'Hit to', 'Overlap length'],
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
	my @features;

	# For all of our features, if we want them, bump coverage #s
	for my $feat_object ($seq_object->get_SeqFeatures) {
		if($wanted_tags{ $feat_object->primary_tag }){
			my $loc = $feat_object->location;
			if(ref($loc) eq 'Bio::Location::Simple'){
				push(@features,[$feat_object->start, $feat_object->end,$bio->_getIdentifier($feat_object),$feat_object]);
			}else{
				my $i = 0;
				for my $location ( $loc->sub_Location ) {
					push(@features,[$location->start, $location->end,$bio->_getIdentifier($feat_object) . '_subpart' . ($i++), $feat_object]);
				}
			}
		}
	}

	# Given that data, calculate overlaps
	foreach my $from_feature ( @features ){
		foreach my $to_feature ( @features ){
			my ($a0,$a1,$ai,$af) = @{$from_feature};
			my ($b0,$b1,$bi,$bf) = @{$to_feature};
			# Some we want to skip so we don't have duplicates
			if($ai eq $bi || $ai ge $bi){
				next;
			}
			if(is_overlapped($a0,$a1,$b0,$b1)){
				push(@data_data, [$ai,$bi, get_overlap_length($a0,$a1,$b0,$b1) ]);
			}
		}
	}
}

sub get_overlap_length {
	my ($a0,$a1,$b0,$b1) = @_;
	my $dist;
	for(my $i=$a0;$i<=$a1;$i++){
		if(is_inside($i,$b0,$b1)){
			$dist++;
		}
	}
	return $dist;
}
sub is_overlapped {
	my ($a0,$a1,$b0,$b1) = @_;
	return
		is_inside($a0,$b0,$b1) ||
		is_inside($a1,$b0,$b1) ||
		is_inside($b0,$a0,$a1) ||
		is_inside($b1,$a0,$a1);
}
sub is_inside {
	my ($idx, $start, $end) = @_;
	return $idx >= $start && $idx <= $end;
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
