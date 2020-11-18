#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Bio::SearchIO;
use Bio::SeqIO;
use File::Temp;
use CPT::GalaxyGetOpt;

# PODNAME: cpt_ppepp.pl

my $ggo = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [ 
		[ 'file|f', 'Input Genbank file', { required => 1, validate => 'File/Input' } ], 
		[ 'lookahead_dist', 'Lookahead distance', { required => 1, default => 30, validate => 'Float', min => '0', max => '100'}],
		[ 'evalue', 'report sequences <= this E-value threshold in output', { validate => 'Float', default => 5}],
		[ 'inc','consider sequences <= this E-value threshold as significant', { validate => 'Float', default => 0.05}],
		[ 'f1', 'Stage 1 (MSV) threshold: promote hits w/ P <= F1', {validate => 'Float', default=> 0.1}],
		[ 'f2', 'Stage 1 (Vit) threshold: promote hits w/ P <= F2', {validate => 'Float', default=> 0.05}],
		[ 'f3', 'Stage 1 (Fwd) threshold: promote hits w/ P <= F3', {validate => 'Float', default=> 0.05}],
		[ 'domz', 'set # of significant seqs, for domain E-value calculation', {validate => 'Float', default=>10}],
	],
	'outputs' => [
		[
			'report',
			'PPEPP Report',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'ppepp',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PPEPP',
		'appname' => 'PPEPP',
		'appdesc' => 'Paradigm Phage Early Promoter Project. Finds promoters that ought to be expressed by unmodified host RNA polymerase, utilising HMMR3 databases curated from T3, T5, T7, P1, P2, P4, and 186',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 't/data/lambda.gbk',
			},
			outputs => { 
				'report' => ["ppepp.Sheet1.csv", "t/data/ppepp.csv"],
			}
		},
	]
);

my %phage_results = ();

use CPT::Analysis::PPEPP;
my $ppepp = CPT::Analysis::PPEPP->new(
	lookahead_dist => $options->{lookahead_dist},
	evalue         => $options->{evalue},
	inc            => $options->{inc},
	f1             => $options->{f1},
	f2             => $options->{f2},
	f3             => $options->{f3},
	domz           => $options->{domz},
);

my @data = $ppepp->analyze_sequence( $options->{file} );

my %results = (
	'Sheet1' => {
		header => ["target name","accession","query name","accession","E-value","score","bias","E-value","score","bias","exp","reg","clu","ov","env","dom","rep","inc","description of target"],
		data => \@data,
	},
);

#use CPT::OutputFiles;
#my $output = CPT::OutputFiles->new(
	#name => 'report',
	#GGO => $ggo,
#);
#$output->CRR(data => \%results);
