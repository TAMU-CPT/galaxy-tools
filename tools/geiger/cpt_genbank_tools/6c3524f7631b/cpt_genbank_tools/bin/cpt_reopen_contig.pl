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

my $libCPT = CPT->new();

my $options = $libCPT->getOptions(
	'options' => [
		[
			'file|f',
			'Input file',
			{
				required => 1,
				validate => 'File/Input',
			}
		],
		[
			'rebase',
'This specified base will be the new first base of the sequence.',
			{ required => 1, validate => 'Int', default => 1000 }
		],
	],
	'outputs' => [
		[
			'genbank',
			'Reopened Genbank File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'reopened',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank',
			}
		],
	],
	'defaults' => [
		'appid'   => 'ReopenContig',
		'appname' => 'Reopen Contig',
		'appvers' => '1.94',
		'appdesc' => 'Reopen Genbank file in a new location',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/split-feature.gbk',
				'rebase' => 3,
			},
			outputs      => {
				'genbank' => ['reopened.gbk', 'test-data/outputs/reopened.gbk'],
			}
		},
	]
);

my $offset        = $options->{'rebase'} - 1;

=head1 LIMITATIONS

Currently limited to single genbank records (no multi-genome genbank files)

=cut

use Bio::SeqIO;
my $seqio_object = Bio::SeqIO->new(-file => $options->{file}, -format=>'genbank');
my $seqobj = $seqio_object->next_seq;
my $genome_length = $seqobj->length;

my $new_sequence =
    $seqobj->subseq( $offset + 1, $genome_length )
  . $seqobj->subseq( 1,           $offset );
$seqobj->seq( $new_sequence, 'DNA' );

# Correct for 1 counting
if ( ($offset) % 3 != 0 ) {
	warn "Offset chosen will cause a reading frame shift";
}

my $moved_feature_count = 0;
my $feature_offset      = 0;

foreach my $feat ( $seqobj->get_SeqFeatures ) {
	my $location_type = ref $feat->location;
	if ( $location_type eq 'Bio::Location::Simple' ) {
		if ( $feat->primary_tag ne 'source' ) {
			move_feature($feat, $genome_length);
		}
	}
	elsif ( $location_type eq 'Bio::Location::Split' ) {
		foreach my $sublocs($feat->location->each_Location()){
			move_location($sublocs, $genome_length);
		}
	}else{
		warn "I don't know how to handle $location_type types of locations. Please submit a bug with Eric Rasche";
	}
}

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'genbank',
	libCPT => $libCPT,
);
$crr_output->CRR(data => $seqobj);


sub move_location {
	my ($loc, $genome_length) = @_;
	my ( $start, $end ) = (
		$loc->start,
		$loc->end
	);
	$start -= $offset;
	$end   -= $offset;
	if ( $start < 1 && $end < 1 ) {
		$start += $genome_length;
		$end   += $genome_length;
	}
	elsif (    ( $start < 1 && $end > 0 )
		|| ( $start > 0 && $end < 1 ) )
	{
		die "Location at which you are attempting to reopen is inside of a feature\nIf you need to be able to disable/allow this, please submit a feature request with Eric";
	}
	$loc->start($start);
	$loc->end($end);
}
sub move_feature {
	my ($feat, $genome_length) = @_;
	my ( $start, $end ) = (
		$feat->location->start,
		$feat->location->end
	);
	$start -= $offset;
	$end   -= $offset;
	if ( $start < 1 && $end < 1 ) {
		$start += $genome_length;
		$end   += $genome_length;
	}
	elsif (    ( $start < 1 && $end > 0 )
		|| ( $start > 0 && $end < 1 ) )
	{
		die "Location at which you are attempting to reopen is inside of a feature\nIf you need to be able to disable/allow this, please submit a feature request with Eric";
	}
	$feat->start($start);
	$feat->end($end);
}


=head1 DESCRIPTION

This tool reopens contigs in a specific location, without otherwise modifying the file. This is B<very important> as the traditional method of doing this has involved adding and removing features in the genbank file. Failure to remove these features when done is B<common>, and results in B<terrible> looking student genome maps. As such, this tool should be used exclusively to all other methods.

=cut
