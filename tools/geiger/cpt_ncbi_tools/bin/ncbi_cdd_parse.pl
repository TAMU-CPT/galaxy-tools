#!/usr/bin/env perl
use strict;
use warnings;

use CPT::GalaxyGetOpt;
my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[
			'file',
			'NCBI CDD Query Output',
			{
				validate => 'File/Input',
			}
		],
	],
	'outputs' => [
		[
			'output',
			'GFF3 formatted results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'out',
				data_format    => 'genomic/interval',
				default_format => 'GFF3',
			}
		],
	],
	'defaults' => [
		'appid'   => 'NCBICDDtoGFF3',
		'appname' => 'NCBI CDD query to GFF3',
		'appdesc' => 'converts results of raw CDD queries to gff3 format for downstream analysis',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name => "Default",
			params => {
				file => 'test-data/inputs/output.txt',
			},
			outputs => {
				'output' => ["out.gff3", 'test-data/outputs/cdd-01.gff3' ],
			},
		},
	],
);

#             0       1          2         3      4    5         6          7           8            9            10            11
my @fields = ('Query','Hit type','PSSM-ID','From','To','E-Value','Bitscore','Accession','Short name','Incomplete','Superfamily','Definition',);

my @features;
open(my $fh, '<', $options->{file});
while(<$fh>){
	chomp;
	next unless ($_ =~ /^Q#/);
	my @row = split(/\t/, $_);
	my $id = $row[0];
	my $desc = "";
	if($row[0] =~ /\[([^ ]*)(.*)\]/){
		$id = $1;
		$desc = $2;
	}
	my $feature = Bio::SeqFeature::Generic->new(
		-seq_id => $id,
		-primary => 'Match',
		-start => $row[3],
		-end => $row[4],
		score => $row[5],
		-frame => '0',
		strand => 0,
		source_tag => 'CDD',
		tag => {
			'Dbxref' => $row[7],
		}
	);
	push(@features, $feature);
}
close($fh);


use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
my $gff_factory = Bio::Tools::GFF->new(-gff_version=>3);
my $gff3 = "##gff-version 3\n";
foreach my $feature(@features){
	$gff3 .= $gff_factory->_gff3_string($feature) . "\n";
}
use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'output',
	GGO => $ggo,
);
$crr_output->CRR(data => $gff3);
