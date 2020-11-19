#!/usr/bin/env perl
#
#       Code written by Eric Rasche
#               mailto:rasche.eric@yandex.ru
#               tel:   404.692.2048
#               http://eric.rasche.co.uk
#       for
#               Center for Phage Technology
#
#
# PODNAME: cpt_analyze-terl.pl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;
use CPT::Analysis::TerL;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		['file'                , 'Input file'                                             , { required => 1           , validate => 'File/Input' } ]    ,
		['blast'               , 'Run Blast analysis']                                    ,
		['hmmer'               , 'Run HMMER analysis']                                    ,
		[]                     ,
		['hmmer_evalue_cutoff' , 'Minimum log evalue acceptable from HMMER'               , {validate => 'Int'        , default => -64                  , required => 1}] ,
		['blast_evalue_cutoff' , 'Minimum log evalue acceptable from BLAST'               , {validate => 'Int'        , default => -140                 , required => 1}] ,
		['blast_dice_cutoff'   , 'Minimum dice score acceptable from BLAST'               , {validate => 'Int'        , default => 30                   , required => 1}] ,
	]                              ,
	'outputs' => [
		[
			'terl_report',
			'Output Results File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'terl_analysis',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid' => 'Analysis_TerL',
		'appname' => 'Analyse TerL Sequences',
		'appdesc' => 'Analyse sequences for similarity to known TerL sequences, assigning a morphology based on homology',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 't/data/test.fa',
				'blast' => '',
				'hmmer' => '',
				'terl_report_format' => 'YAML',
			},
			outputs      => {
				'terl_report' => ['terl_analysis.yml','t/data/output.yml'],
			}
		},
	],
);

my $terl = CPT::Analysis::TerL->new(
	hmmer_evalue_cutoff => $options->{'hmmer_evalue_cutoff'},
	blast_evalue_cutoff => $options->{'blast_evalue_cutoff'},
	blast_dice_cutoff   => $options->{'blast_dice_cutoff'},
	search_hmmer        => $options->{'hmmer'},
	search_blast        => $options->{'blast'},
);

my $output_data_ref = $terl->run(
	input_file => $options->file,
);

use CPT::OutputFiles;
my $out = CPT::OutputFiles->new(
	name => 'terl_report',
	GGO => $ggo,
);
$out->CRR(data => $output_data_ref);
