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
use Storable qw/dclone/;

my $libCPT  = CPT->new();
my $options = $libCPT->getOptions(
	'options' => [
		[ 'file|f', 'Input file',
			{
				validate => 'File/Input',
				required => 1,
			}
		],
		[ 'end', 'End of Terminal Repeat Region. READ THE DIRECTIONS BELOW (or in POD)!!!!!!!',
			{
				validate => 'Int',
			}
		],
	],
	'outputs' => [
		[
			'results',
			'Genbank File with duplicated TRs',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'with_terminal_repeats',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank',
			}
		],
	],
	'defaults' => [
		'appid'   => 'Add_TRs',
		'appname' => 'Add Terminal Repeats',
		'appdesc' => 'to a Genbank file with one or more genomes',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Repeat Region",
			params => {
				'file' => 'test-data/inputs/repeat_region.gbk',
			},
			outputs => {
				'results' => [ 'with_terminal_repeats.gbk', 'test-data/outputs/add_tr.repeat_region_at_0.gbk' ]
			},
		},
		{
			test_name    => "Repeat Region Offset",
			params => {
				'file' => 'test-data/inputs/repeat_region2.gbk',
			},
			outputs => {
				'results' => [ 'with_terminal_repeats.gbk', 'test-data/outputs/add_tr.repeat_region_not_0.gbk' ]
			},
		},
	],
);

use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;
my $seqio_object = Bio::SeqIO->new(-file => $options->{file}, -format=>'Genbank');
my $number_of_genomes = 0;
my $seq_object = $seqio_object->next_seq;

my $repeat_start = 1;
my $repeat_end = 1;
my $repeats_encountered = 0;
my $repeat_feat;
if(! defined($options->{end})){
	foreach my $feat($seq_object->get_SeqFeatures){
		if($feat->primary_tag eq 'repeat_region'){
			$repeats_encountered++;
			if($repeats_encountered > 1){
				die 'Script encountered more than one repeat_region feature, remove all others first!';
			}
			$repeat_start = $feat->start();
			$repeat_end = $feat->end();
			$repeat_feat = dclone($feat);
		}
	}
}else{
	$repeat_end = $options->{end};
	$repeats_encountered++;
}

print STDERR "Cutting from $repeat_start to $repeat_end\n";
# Grab region from 1 .. end
my $tr = $seq_object->subseq($repeat_start, $repeat_end);
my $len = $seq_object->length(); # this is fine since features + TR region are just appended to the end of the genome
print STDERR "Length is now $len\n";
my $new_seq = $seq_object->seq() . $tr;
# Add our TR to it

$seq_object->seq($new_seq);

my @clonefeats = map { dclone($_) } get_features_in_region($seq_object, $repeat_start, $repeat_end);

for my $clonefeat(@clonefeats, $repeat_feat){
	$clonefeat->display_name( $clonefeat->display_name() . '_rep');
	if($clonefeat->has_tag('locus_tag')){
		my @loci = map { $_ . '_rep' } $clonefeat->get_tag_values('locus_tag');
		$clonefeat->remove_tag('locus_tag');
		$clonefeat->add_tag_value('locus_tag', @loci);
	}
	$clonefeat->start($clonefeat->start() + $len - $repeat_start + 1);
	$clonefeat->end($clonefeat->end() + $len - $repeat_start + 1);
	$seq_object->add_SeqFeature($clonefeat);
}

sub get_features_in_region {
	my ($seq_object, $start, $end) = @_;
	my @feats;
	foreach my $feat($seq_object->get_SeqFeatures){
		if($feat->start >= $repeat_start && $feat->end < $repeat_end){
			push(@feats, $feat);
		}
	}
	return @feats;
}


use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
	name => 'results',
	libCPT => $libCPT,
);
$output->CRR(data => $seq_object);


=head1 NAME

TR Addition

=head1 DESCRIPTION

This tools adds TR regions to a Genbank file.

The region you supply duplicated and added to the end of the genome. The features from that region are likewise duplicated, and "_rep" added to their locus tag to indicate repeated features. If the specified TR region ends inside a feature, it is not duplicated. This tool does B<not> renumber.

=head1 REQUIREMENTS

Your genomes MUST be opened such that the start of the terminal repeat is the first base! This is unfortunately non-negotiable. If the repeat_region is in the middle of the genome, it will be duplicated as-is, and appended to the end. This feature can be abused, but you B<must> be aware of it. See the section L</"Feature Abuse">

You may specify repeats in one of the following ways:

=over 4

=item Add a repeat_region feature to the start of your genome in artemis

=item Specify a repeat end location

=back

Either of these methods is sufficient to add terminal repeats.


=head1 Feature Abuse

Sometimes you may find yourself in the following situation:

=over 4

=item You have a genome with a C<repeat_region>

=item During sequencing, this repeat region was not fully collapsed by the sequencer

=item The full repeat is at the start

=item A portion of the front of the repeat is also at the end (maybe this was done to improve annotation?)

=back

In this event, you need to do the following:

=over 4

=item Calculate how much of the repeat region was duplicated

=item Create the C<repeat_region> feature at the start of the genome

=item Move the feature up by the number of duplicated bases

=back

Then, when the tool is run, this subset of the (true) repeat region will be duplicated and features added to the end. After that, you will need to fix the location of both repeat regions to reflect the true repeat region

=cut
