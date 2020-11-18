#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use CPT::Bio::ORF;
use CPT::Bio::SAR;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file', 'Input file',
			{ required => 1, validate => 'File/Input' }
		],
		[
			'min_gene',
			'Minimum gene length for ORF finding step',
			{ required => 1, validate => 'Int', default => 30, }
		],
	],
	'outputs' => [
		[
			'fasta',
			'Found ORFs',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'orfs',
				data_format    => 'genomic/raw',
				default_format => 'Fasta',
			}
		],
	],
	'defaults' => [
		'appid'   => 'SarFinder',
		'appname' => 'Sar Finder',
		'appvers' => '1.96',
		'appdesc' => 'Attempts to locate SAR domains according to sequence rules',
	],
	'tests' => [
	],
);


use CPT::Bio::ORF;
my $orf = CPT::Bio::ORF->new(
	min_gene_length => $options->{min_gene},
);
use CPT::Bio;
my $bio = CPT::Bio->new();
my $seq_obj = ${ $bio->requestCopy( 'file' => $options->{file} ) };
my @seqs = $orf->run($seq_obj->seq());

my $sar = CPT::SAR->new();
my @sar_seqs = $sar->filter_sar(@seqs);

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'fasta',
	GGO => $ggo,
);
$crr_output->CRR(data => \@sar_seqs);

=head1 DESCRIPTION

Runs the CPT's SAR finder to locate possible SAR domains

=cut
