#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my %options = ( a => "Option A", b => "Option B" );

my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[
			'file', 'Input file',
			{
				validate => 'File/Input',
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
		[
			"option" => "Select an option!",
			{
				validate => 'Option',
				options  => \%options,
				multiple => 1,
			}
		],
		[
			"float" => "I'm a float",
			{ validate => 'Float' }
		],
		[
			"int" => "I'm an int",
			{ validate => 'Int', default => [42, 34], required => 1, multiple => 1 }
		],
		[],
		['New section'],
		[
			"string" => "I'm a simple string",
			{ validate => 'String' }
		],
		[
			'flag' => 'This is a flag',
			{ validate => 'Flag' }
		],
	],
	'outputs' => [
		[
			'my_output_data',
			'Output TXT File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'out',
				data_format    => 'text/plain',
				default_format => 'TXT',
			}
		],
	],
	'defaults' => [
		'appid'   => 'TemplateScript',
		'appname' => 'Template',
		'appdesc' => 'A basic template for a new CPT application',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
			},
			outputs => {
				'my_output_data' => ["out.txt", 'test-data/outputs/template.default.txt' ],
			},
		},
		{
			test_name    => "Option A specified",
			params => {
				'int', '10000',
			},
			outputs => {
				'my_output_data' => ["out.txt", 'test-data/outputs/template.A.txt' ],
			},
		},
	],
);

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'my_output_data',
	GGO => $ggo,
);
$crr_output->CRR(data => sprintf( "Hello, World\n%s\n", join(', ', @{$options->{int}}) ));

=head1 NAME

TemplateScript

=head1 DESCRIPTION

TemplateScript does all sorts of fun things like twiddling twaddles and widgeting wobs. Now it even bobs flooby bits!

=cut
