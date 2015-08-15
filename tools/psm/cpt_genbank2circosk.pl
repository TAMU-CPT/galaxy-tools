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
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();
my %colors = map { $_ => $_ } qw(red orange yellow green blue gray black white);
my %intensity = map { $_ => $_ } qw (vvvvl vvvl vvl vl vd vvd vvvd vvvvd);

my $options = $ggo->getOptions(
	'options' => [
		[
			'file',
			'Input file',
			{
				required => 1,
				validate => 'File/Input',
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
		[
			'chromosome' => 'Name for the
			chromosome inside Circos',
			{
				required => 1,
				validate => 'String'
			}
		],
		[
			'color' => 'Color to use for Circos plot',
			{
				required => 1,
				validate => 'Option',
				options  => \%colors,
			}
		],
		[
			'intensity' => 'Circos color intensity. ',
			{
				validate => 'Option',
				options  => \%intensity,
			}
		],
	],
	'outputs' => [
		[
			'circosk',
			'Circos Karyotype File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'karyotype',
				data_format    => 'text/plain',
				default_format => 'TXT'
			}
		],
	],
	'defaults' => [
		'appname' => 'Genbank2CircosK',
		'appid'   => 'Genbank2CircosK',
		'appvers' => '1.94',
		'appdesc' =>
'Convert genbank files to Circos Karyotype configuration files',
	],
	'tests' => [
		{
			test_name => "Default",
			params => {
				'file' => 'test-data/inputs/multi.gbk',
				'chromosome' => 'test',
				'color' => 'red',
				'intensity' => 'vvvl',
			},
			outputs => {
				'circosk' => ['karyotype.txt', 'test-data/outputs/circosk.conf'],
			}
		},
	],
);

use CPT::Bio;
my $bio = CPT::Bio->new();

my @results;
my $c = 0;
my $seqio = $bio->requestCopyIO( file => $options->{file} );

while(my $seqobj = $seqio->next_seq()){
	foreach my $feat ( $seqobj->get_SeqFeatures () ) {

		#band test 12 CDS__test_1gbk 5715 6335 red]
		next if ( $feat->primary_tag ne 'CDS' );
		my $id = $bio->_getIdentifier($feat);
		$id =~ s/\s+/_/g;
		push(
			@results,
			join(
				' ',
				(
					'band',
					$options->{'chromosome'},
					$c++,
					$id,
					$feat->start,
					$feat->end,
					(defined $options->{'intensity'} ? $options->{'intensity'} : '') . $options->{'color'}
				)
			)
		);
	}
}

my $z = join( "\n", @results );

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
        name => 'circosk',
        GGO => $ggo,
);
$output->CRR(data => $z);
