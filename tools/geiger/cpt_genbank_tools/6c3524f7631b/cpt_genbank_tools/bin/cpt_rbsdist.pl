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
			'RBS Distance Results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'rbsdist',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'RBSDist',
		'appname' => 'RBS Distance',
		'appdesc' => 'distance between CDS and associated RBS feature',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);

# Which tags should we consider as part of PIGS? These should be coding sequences
my %wanted_tags = (
	CDS => 1,
);

my %data = (
	'Sheet1' => {
		header => ['Sequence', 'Feature', 'Associated RBS', 'RBS Start', 'RBS End', 'Distance', 'RBS Sequence'],
		data => [],
	}
);
my @data_data;

use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;
my $seqio_object = Bio::SeqIO->new(-file => $options->{file}, -format=>'genbank');
# For all genomes in the GBK file
while(my $seq_object = $seqio_object->next_seq){
	my $seqid = $seq_object->display_id();
	my @rbs_map;
	foreach my $rbs($seq_object->get_SeqFeatures()){
		if($rbs->primary_tag eq 'RBS'){
			push(@rbs_map, $rbs);
		}
	}
	foreach my $feat($seq_object->get_SeqFeatures()){
		if($feat->primary_tag eq 'CDS'){
			my $rbs = locate_upstream_rbs($feat,\@rbs_map);
			if(defined $rbs){
				push(@data_data, [
					$seqid,
					$bio->_getIdentifier($feat),
					$bio->_getIdentifier($rbs),
					$rbs->start,
					$rbs->end,
					inter_feature_spread($feat, $rbs),
					$bio->get_seq_from_feature($rbs),
				]);
			}else{
				push(@data_data, [
					$seqid,
					$bio->_getIdentifier($feat),
					'NONE FOUND',
					''
				]);
			}
		}
	}
}


sub locate_upstream_rbs {
	my ($feat, $rbs_ref) = @_;
	my @rbs_list = @{$rbs_ref};
	my @close;
	foreach my $rbs(@rbs_list){
		if(is_close($feat, $rbs)){
			push(@close, $rbs);
		}
	}
	if(scalar @close > 1){
		warn "Found MORE THAN ONE possible RBS for " . $bio->_getIdentifier($feat) . ": " . join("\n", map { "\t" . $bio->_getIdentifier($_) } @close ) . "\n";
	}
	return $close[0];
}

sub is_close{
	my ($feat_a, $feat_b) = @_;
	if($feat_a->strand() != $feat_b->strand()){
		return 0;
	}
	if(inter_feature_spread($feat_a, $feat_b) < 15){
		return 1;
	}
	return 0;
}

sub inter_feature_spread {
	my ($feat_a, $feat_b) = @_;
	if($feat_a->strand() != $feat_b->strand()){
		warn "Cannot compare distance for features on different strands";
		return -1;
	}
	if($feat_a->strand() == 1){
		return abs($feat_a->start() - $feat_b->end()) - 1;
	}
	if($feat_a->strand() == -1){
		return abs($feat_a->end() - $feat_b->start()) - 1;
	}
	return abs($feat_a->start() - $feat_b->end());
}






$data{Sheet1}{data} = \@data_data;
use CPT::OutputFiles;
my $csv_output = CPT::OutputFiles->new(
        name => 'results',
        libCPT => $libCPT,
);
$csv_output->CRR(data => \%data);

=head1 DESCRIPTION

RBS Distance calculator. Produces a table of the separation, location, and sequence of the upstream RBS from every feature.

=cut
