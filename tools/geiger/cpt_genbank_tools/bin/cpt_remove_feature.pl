#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file|f',
			'Input file',
			{
				required => 1,
				validate => 'File/Input',
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
		[
			'feature',
			'Name of the feature to be removed',
			{
				required => 1,
				multiple => 1,
				validate => 'Genomic/Tag'
			}
		],
	],
	'outputs' => [
		[
			'genbank',
			'Removed Tag',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'removed_tag',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank',
			}
		],
	],
	'defaults' => [
		'appid'   => 'RemoveFeature',
		'appname' => 'Remove Feature',
		'appvers' => '1.94',
		'appdesc' => 'Remove a feature from a Genbank file',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				"file" => "test-data/inputs/multi.gbk",
				"feature" => 'CDS',
			},
			outputs => {
				"results" => ["removed_tag.gbk", "test-data/outputs/removed_cds.gbk"]
			},
		},
	],
);

my %args     = ( 'file' => $options->file, );

use CPT::Bio;
my $bio = CPT::Bio->new();
my $seqobj   = ${ $bio->requestCopy(%args) };
my %bad_keys = map { $_ => 1 } @{ $options->{'feature'} };

my $saved   = 0;
my $removed = 0;
my @good_features;
foreach my $feat ( $seqobj->get_SeqFeatures ) {
	if ( !$bad_keys{ $feat->primary_tag } ) {
		$saved++;
		push( @good_features, $feat );
	}
	else {
		$removed++;
	}
}

$seqobj->remove_SeqFeatures();
foreach (@good_features) {
	$seqobj->add_SeqFeature($_);
}

#print "Saved $saved features. Removed $removed features\n";
use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'genbank',
	GGO => $ggo,
);
$crr_output->CRR(data => $seqobj);

=head1 DESCRIPTION

Batch removal of features from a GenBank file

=cut