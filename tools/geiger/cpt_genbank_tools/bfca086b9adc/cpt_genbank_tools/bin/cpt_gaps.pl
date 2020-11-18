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
		'appdesc' => 'calculates statistics on gaps in genomic sequences',
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
		header => ['Sequence', 'Sum', 'Mean', 'SD', 'Min', 'Max', 'Median' ,'Corrected Mean', 'Corrected Median'],
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
		if($feat_object->primary_tag ne 'source'){
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

	use Statistics::Descriptive;
	my $stat = Statistics::Descriptive::Full->new();

	my $last_state = $coverage{1};
	my $string='';
	for(my $i = 1; $i<$seq_object->length();$i++){
		$string .= $coverage{$i};
	}

	while($string =~ /([0]+)/g){
		$stat->add_data(abs($+[0] - $-[0]));
	}
	push(@data_data, [$seq_object->display_id(), $stat->sum(), $stat->mean(), $stat->standard_deviation(), $stat->min(), $stat->max(), $stat->median(), 100*$stat->mean()/$seq_object->length(), 100*$stat->median()/$seq_object->length()]);
}

$data{Sheet1}{data} = \@data_data;
use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'results',
	libCPT => $libCPT,
);
$crr_output->CRR(data => \%data);

=head1 NAME

Genomic OverlAp caLculator (GOAL)

=head1 DESCRIPTION

Locate regions of overlap between genes.

=cut
