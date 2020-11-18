#!/usr/bin/perl
#
#       Code written by Eric Rasche
#               mailto:rasche.eric@yandex.ru
#               tel:   404.692.2048
#               http://eric.rasche.co.uk
#       for
#               Center for Phage Technology
#

use strict;
use warnings;

use CPT;
use Data::Dumper;

my $libCPT  = CPT->new();
my $options = $libCPT->getOptions(
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
			my $gapstr = cigar_from_alns($hsp->query_string, $hsp->homology_string, $hsp->hit_string);
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
					'Gap' => $gapstr,
				},
				
			);
			push(@features, $feature);
		}
	}
}

sub cigar_from_alns {
	my ($q, $m, $h) = @_;
	my $strict_m = 1;

	use List::MoreUtils qw/each_array/;
	my @q = split //,$q;
	my @m = split //,$m;
	my @h = split //,$h;
	my $it = each_array(@q,@m,@h);

	my $last = "None";
	my $last_count = 0;
	my $cigar_line;
	while(my ($qi,$mi,$si) = $it->()){
		my $chr = '';
		# Match or similar
		if($mi ne ' ' || $mi eq '+'){
			# d = skip in q
			$chr = '=';
		}
		elsif($mi eq ' '){
			if($qi eq '-'){
				$chr = 'D';
			}elsif($si eq '-'){
				$chr = "I";
			}else{
				$chr = 'X';
			}
		}
		else{
			warn 'bad data';
		}

		if($strict_m){
			if($chr eq '=' || $chr eq 'X'){
				$chr = 'M';
			}
		}

		$last_count += 1;
		
		if($chr ne $last){
			if($last ne 'None'){
				$cigar_line .= "$last$last_count ";
			}
			$last_count = 0;
			$last = $chr;
		}
	}
	if($last ne 'None'){
		$cigar_line .= "$last$last_count ";
	}
	return $cigar_line;
}

my $gff3 = "##gff-version 3\n";
foreach my $feature(@features){
	$gff3 .= $gff_factory->_gff3_string($feature) . "\n";
}
use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'output',
	libCPT => $libCPT,
);
$crr_output->CRR(data => $gff3);

=head1 DESCRIPTION

BlastXML files, when transformed to GFF3, do not normally show gaps in the blast hits. This tool aims to fill that "gap".

For an input BlastXML file, this tool will produce a GFF3 file containing all of the relevant information: start, stop, score, and the GAP as a CIGAR string.

L<One CIGAR Spec|http://www.sequenceontology.org/gff3.shtml> does not list X and = for mismatches and matches, respsectively. However, L<Other CIGAR Specs|http://samtools.github.io/hts-specs/SAMv1.pdf> do allow for those characters. Given that I cannot anticipate which characters your downstream analysis methods will support, this tool provides the option to use just M (C<strict_m> enabled), or to use both, as is the default.

=cut
