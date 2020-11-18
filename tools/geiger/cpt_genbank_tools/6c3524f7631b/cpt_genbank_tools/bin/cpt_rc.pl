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
	],
	'outputs' => [
		[
			'results',
			'RC\'d Genbank File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'rc',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank',
			}
		],
	],
	'defaults' => [
		'appid'   => 'RevCom',
		'appname' => 'Reverse and Complement',
		'appdesc' => 'for a given genbank file',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default GBK",
			params => {
				'file' => 'test-data/inputs/single.gbk',
			},
			outputs => {
				results => ['rc.gbk', 'test-data/outputs/rc.gbk'],
			}
		},
	],
);

use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;
use Bio::Seq;
my $seqio_object = Bio::SeqIO->new(-file => $options->{file}, -format=>'Genbank');
my $number_of_genomes = 0;
my @new_seqs;
while(my $seq_object = $seqio_object->next_seq){
	my $seq = $seq_object->seq();
	my $len = $seq_object->length();
	# We're going to be creating a new sequence object pretty much no
	# matter what. might as well be a BP one with the revcom method.
	my $bs = Bio::Seq->new( -seq => $seq, -display_id => "asdf");
	$seq_object->seq($bs->revcom->seq());

	# For all of our features, if we want them, bump coverage #s
	for my $feat_object ($seq_object->get_SeqFeatures) {
		my $loc = $feat_object->location;
		if(ref($loc) eq 'Bio::Location::Simple'){
			$loc->start($len - $loc->start() + 1);
			$loc->end($len - $loc->end() + 1);
		}else{
			for my $location ( $loc->sub_Location ) {
				$location->start($len - $location->start() + 1);
				$location->end($len - $location->end() + 1);
			}
		}
		
		$feat_object->strand(-$feat_object->strand());
	}
	push(@new_seqs, $seq_object);
}

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
	name => 'results',
	libCPT => $libCPT,
);
$output->CRR(data => \@new_seqs);

=head1 NAME

RevCom/Reverse and Complement

=head1 DESCRIPTION

Simply flips the sequence (+comlements), and does the same for features

=cut
