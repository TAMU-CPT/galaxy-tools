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
				validate => 'File/Input'
			}
		],
	],
	'outputs' => [
		[
			'report',
			'List of Hidden Stops ',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'hidden',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'HiddenStopDetector',
		'appname' => 'Hidden Stop Codon Detector',
		'appvers' => '1.94',
		'appdesc' => 'detects +1/-1 stop codons',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/single.gbk',
			},
			outputs      => {
				'report' => ['hidden.Sheet1.csv', 'test-data/outputs/hidden.csv'],
			}
		},
	],
);

use CPT::Bio;
my $bio = CPT::Bio->new();

my %args = (
	'file'     => $options->file,
	'callback' => \&func,
	'header'   => 1,
	'subset'   => 'CDS',
);
$bio->parseFile(%args);


sub func {
	my $response_ref = shift;
	my @response     = @$response_ref;
	my $answer;

	my @data;
	my %stop_codons = map { $_ => 1 } qw(TAG TAA TGA);
	foreach (@response) {
		my ( $header, $seq ) = @{$_};
		for ( my $i = 0 ; $i < length($seq) - 2 ; $i += 3 ) {

			# +1
			if ( $stop_codons{ substr( $seq, $i + 1, 3 ) } ) {
				push( @data, [ $header, '+1', $i + 2 ] );
			}
			if ( $stop_codons{ substr( $seq, $i + 2, 3 ) } ) {
				push( @data, [ $header, '-1', $i + 3 ] );
			}

			# +2
		}
	}

	my %results = (
		'Sheet1' => {
			'header' => [ "Parent Sequence", "Frame", "Location" ],
			'data'   => \@data,
		},
	);

	use CPT::OutputFiles;
	my $crr_output = CPT::OutputFiles->new(
		name => 'report',
		libCPT => $libCPT,
	);
	$crr_output->CRR(data => \%results);
}

=head1 DESCRIPTION

L<Hidden Stops|https://en.wikipedia.org/wiki/Stop_codon#Hidden_stops> "are non-stop codons that would be read as stop codons if they were frameshifted +1 or -1."

=cut
