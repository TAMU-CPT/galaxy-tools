#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Bio::Tools::SeqStats;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file' => 'Input file',
			{
				validate => 'File/Input',
				required => 1,
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
	],
	'outputs' => [
		[
			'converted',
			'Converted Output File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'converted',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank',
			}
		],
	],
	'defaults' => [
		'appid'   => 'GenomicConverter',
		'appname' => 'File Format Convert',
		'appvers' => '1.94',
		'appdesc' => 'Convert amongst support BioPerl/libCPT formats',
	],
	'tests' => [
	],
);
use CPT::Bio;
my $bio = CPT::Bio->new();

use Bio::SeqIO;
my $seqio_object = $bio->getSeqIO($options->{file});

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'converted',
	GGO => $ggo,
);
$crr_output->CRR(data => $seqio_object);

=head1 DESCRIPTION

This tool converts between file formats supported by BOTH BioPerl AND LibCPT.

=head2 Input Formats

Right now this consists of Genbank and EMBL. If you need more, create a gitlab bug for them

=head2 Output Formats

Most of the L<BioPerl formats|http://bioperl.org/wiki/HOWTO:SeqIO> are supported. The support of C<scf, abi, alf, pln, exp, ctf, ztr> formats are dependent on the server having C<bioperl-ext> installed.

=cut
