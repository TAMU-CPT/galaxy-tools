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

my $libCPT = CPT->new();

my $options = $libCPT->getOptions(
	'options' => [
		[
			'file_mga' => 'Input MGA formatted calls',
			{
				required    => 1,
				validate    => 'File/Input',
			}
		],
	],
	'outputs' => [
		[
			'gff3',
			'GFF3 Formatted Output',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'mga',
				data_format    => 'genomic/interval',
				default_format => 'GFF3'
			}
		],
	],
	'defaults' => [
		'appid'   => 'MGA.convertoGFF3',
		'appname' => 'MGA to GFF3',
		'appvers' => '1.94',
		'appdesc' =>
		  'converts the output of MetaGeneAnnotator to GFF3',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/mga.txt',
			},
			outputs      => {
				'results' => ['mga.gff3', 'test-data/outputs/mga.gff3'],
			}
		},
	]
);

# Load the fasta
use Bio::SeqIO;
use CPT::Bio;
my $bio = CPT::Bio->new();
my $fasta_seqio = $bio->getSeqIO($fasta);
my ( $seq_id, $seq );
while ( my $seqobj = $fasta_seqio->next_seq() ) {
	$seq_id = $seqobj->display_id();
	$seq    = $seqobj->seq();
}

# New seqio;
use Bio::SeqFeature::Generic;
my @features;

# Load/open GFF crap
open( my $fh, '<', $options->{file_mga} );
# Loop over features
while ( <$fh> ){
	next if ($_ =~ /^#/);
	chomp $_;
	my ($gene_id, $start, $end, $strand, $phase, $complete, $score, $model, $rbs_start, $rbs_end, $rbs_score) = split(/\t/, $_);

	$gene_id =~ s/gene_//g;

	my ($rbs, $cds, $gene);

	if($rbs_start ne '-'){
		$rbs = new Bio::SeqFeature::Generic(
			-start  => $rbs_start,
			-end    => $rbs_end,
			-strand => ($strand eq '+' ? 1 : -1),
			-score  => $rbs_score,
			-primary_tag => 'RBS',
			-tag    => {
				ID        => "rbs$gene_id",
				Parent    => "gene$gene_id",
				locus_tag => "gene$gene_id",
			}
		);
	}
	$cds = new Bio::SeqFeature::Generic(
		-start  => $start,
		-phase  => $phase,
		-end    => $end,
		-strand => ($strand eq '+' ? 1 : -1),
		-score  => $score,
		-primary_tag => 'CDS',
		-tag    => {
			ID        => "cds$gene_id",
			Parent    => "gene$gene_id",
			locus_tag => "gene$gene_id",
		}
	);

	$gene = new Bio::SeqFeature::Generic(
		-start  => ($rbs_start ne '-' ? $rbs_start : $start + $phase),
		-end    => $end,
		-strand => ($strand eq '+' ? 1 : -1),
		-primary_tag => 'gene',
		-tag    => {
			ID        => "gene$gene_id",
			locus_tag => "gene$gene_id",
		}
	);
	push(@features, $gene, $cds);
	if(defined $rbs){
		push(@features, $rbs);
	}
}

use Bio::Tools::GFF;
my $gff_factory = Bio::Tools::GFF->new(-gff_version=>3);
my $gff3 = "##gff-version 3\n";
foreach my $feature(@features){
	$gff3 .= $gff_factory->_gff3_string($feature) . "\n";
}

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
	name => 'gff3',
	libCPT => $libCPT,
);
$output->CRR(data => $gff3);
