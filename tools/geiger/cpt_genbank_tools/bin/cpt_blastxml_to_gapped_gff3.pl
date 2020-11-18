#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[
			'file',
			'Input file',
			{
				validate => 'File/Input',
				file_format => ['blastxml'],
			}
		],
		[ 'strict_m', 'One CIGAR spec specificies that Matches AND Mismatches are both represented as M. Another spec allows for = and X to disambiguate. Depending on your downstream pipeline, choose to use =/X or not.', { validate => 'Flag' }],
	],
	'outputs' => [
		[
			'output',
			'Output GFF3 Gapped Alignment',
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
		'appid'   => 'BlastXMLtoGFF3',
		'appname' => 'BlastXML to GFF3 Gapped Alignment',
		'appdesc' => 'converts blast XML data to gapped alignment',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name => "Default",
			params => {
				file => 'test-data/inputs/ex.blastxml',
			},
			outputs => {
				'output' => ["out.gff3", 'test-data/outputs/blast2gff.gff3' ],
			},
		},
	],
);

use Bio::Tools::GFF;
use Bio::SeqFeature::Generic;
my $gff_factory = Bio::Tools::GFF->new(-gff_version=>3);
my @features;

use Bio::SearchIO;
my $in = new Bio::SearchIO(-format => 'blastxml', -file => $options->{file});
my $match_idx = 0;
while( my $result = $in->next_result ) {
	## $result is a Bio::Search::Result::ResultI compliant object
	while( my $hit = $result->next_hit ) {
		## $hit is a Bio::Search::Hit::HitI compliant object
		while( my $hsp = $hit->next_hsp ) {
			my $feature = Bio::SeqFeature::Generic->new(
				-seq_id  => $result->query_name,
				-primary => 'Match',
				-start   => $hsp->start,
				-end     => $hsp->end,
				-score   => $hsp->evalue,
				-frame => '0',
				-source_tag => $hsp->source_tag,
				# taken from bp_search2gff
				-strand      => $hit->strand('query') * $hit->strand('hit'),
				-tag     => {
					'ID' => sprintf('match%05d',++$match_idx),
					'Gap' => $hsp->cigar_string(),
				},
				
			);
			push(@features, $feature);
		}
	}
}

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

=head1 DESCRIPTION

BlastXML files, when transformed to GFF3, do not normally show gaps in the blast hits. This tool aims to fill that "gap".

For an input BlastXML file, this tool will produce a GFF3 file containing all of the relevant information: start, stop, score, and the GAP as a CIGAR string.

L<One CIGAR Spec|http://www.sequenceontology.org/gff3.shtml> does not list X and = for mismatches and matches, respsectively. However, L<Other CIGAR Specs|http://samtools.github.io/hts-specs/SAMv1.pdf> do allow for those characters. Given that I cannot anticipate which characters your downstream analysis methods will support, this tool provides the option to use just M (C<strict_m> enabled), or to use both, as is the default.

=cut
