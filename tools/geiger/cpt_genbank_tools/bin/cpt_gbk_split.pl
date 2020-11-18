#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file',
			'File to split into component genbank files',
			{ required => 1, validate => 'File/Input', 
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
	],
	'outputs' => [
		[
			'genbank',
			'Output Genbank Files',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'downloaded',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank'
			}
		],
	],
	'defaults' => [
		'appid'   => 'cpt.gbk.split',
		'appname' => 'Genbank Split',
		'appdesc' => 'splits a multi-genome genbank file into multi gbk files',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/multi.gbk'
			},
			outputs => {
				'unknown_1.gbk' => ["unknown_1.gbk", 'test-data/outputs/unknown_1.gbk' ],
				'unknown_2.gbk' => ["unknown_2.gbk", 'test-data/outputs/unknown_2.gbk' ],
			},
		},
	]
);


use CPT::OutputFiles;
my $gbk_output = CPT::OutputFiles->new(
	name => 'genbank',
	GGO => $ggo,
);

use CPT::Bio;
my $bio = CPT::Bio->new();

my $idx = 0;
foreach my $file ($options->{file}){
	my $seqio = $bio->getSeqIO($options->{file});
	while(my $seqobj = $seqio->next_seq){
		# Wow this is simple.
		my $name = $seqobj->accession_number();
		$name =~ s/_//g;
		$gbk_output->varCRR(data => $seqobj, filename => sprintf("%s_%s", $name, ++$idx));
	}
}
