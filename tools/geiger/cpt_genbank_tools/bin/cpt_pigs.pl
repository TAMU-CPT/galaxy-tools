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
			'PIGS Results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'pigs',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PIGS',
		'appname' => 'PIGS',
		'appdesc' => 'calculates sum length of intergenic regions as a percentage of the genome',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/multi.gbk',
			},
			outputs => {
				results => ['pigs.Sheet1.csv', 'test-data/outputs/pigs.csv'],
			}
		},
	],
);

# Which tags should we consider as part of PIGS? These should be coding sequences
my %wanted_tags = (
	CDS => 1,
);

my %data = (
	'Sheet1' => {
		header => ['Sequence', 'Coverage (Artemis)', 'Coverage (PIGS)'],
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
	# Coverage history, we'll count 0s vs 1s after
	my @coverage;

	# For all of our features, if we want them, bump coverage #s
	for my $feat_object ($seq_object->get_SeqFeatures) {
		if($wanted_tags{ $feat_object->primary_tag }){
			my $loc = $feat_object->location;
			if(ref($loc) eq 'Bio::Location::Simple'){
				for(my $i=$feat_object->start; $i < $feat_object->end; $i++){
					$coverage[$i]++;
				}
			}elsif(ref($loc) ne 'Bio::Location::Fuzzy'){
				for my $location ( $loc->sub_Location ) {
					for(my $i=$location->start; $i < $location->end; $i++){
						$coverage[$i]++;
					}
				}
			}
		}
	}

	# Given that data, calculate coverages.
	my ($start, $end) = (1, $seq_object->length());
	my ($pigs_covered, $art_covered) = (0,0);
	for(my $i=$start;$i<$end;$i++){
		if(defined $coverage[$i] && $coverage[$i] > 0){
			$pigs_covered++;
			$art_covered += $coverage[$i];
		}
	}
	# Copy data to 2D array
	push(@data_data, [ $seq_object->display_id(), 100*$art_covered/$end, 100*$pigs_covered/$end]);
}

$data{Sheet1}{data} = \@data_data;

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
        name => 'results',
        GGO => $ggo,
);
$output->CRR(data => \%data);

=head1 NAME

PIGS

=head1 DESCRIPTION

Determines percent of intergenic sequence relative to the whole genome. This script was intended to complement the other method for calculating coding density.

There are two methods of calculating coding density:

=over 4

=item C<( sum of lengths of all coding sequences ) / ( length of genome )>

=item C<( sum of lengths of regions without coding sequences ) / ( length of genome )>

=back

the CPT favors the second method as it more accurately reflects the portion of the genome that is not covered in genes. This is relevant as phages have very high coding densities, and a high PIGS score indicates there may be a region which is missing a gene.


=cut
