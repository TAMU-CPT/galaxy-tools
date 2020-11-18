#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Bio::Seq;
use Bio::SeqIO;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file|f',
			'Input file',
			{
				required => 1,
				validate => 'File/Input',
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
		[
			"tag",
			"Analyse by a specific tag",
			{ validate => 'Genomic/Tag', default => 'CDS' }
		],
		[
			'lookahead_max',
			'Maximum separation between RBS and first codon',
			{
				required => 1,
				validate => 'Int',
				default  => 12
			}
		],
		[
			'lookahead_min',
			'Minimum separation between RBS and first codon',
			{
				required => 1,
				validate => 'Int',
				default  => 5
			}
		],
	],
	'outputs' => [
		[
			'report',
			'Shinefind Report',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'shinefind',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'Shinefind',
		'appname' => 'Shinefind',
		'appvers' => '1.99',
		'appdesc' =>
'Create a table of shine-delgarno sequence upstream of proteins',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/single.gbk',
				'lookahead_max' => 30,
				'lookahead_min' => 0,
			},
			outputs      => {
				'report' => ["shinefind.Sheet1.csv", 'test-data/outputs/shinefind.csv'],
			}
		},
	],
);

use CPT::Bio;
my $bio     = CPT::Bio->new();
my %results = (
	'Sheet1' => {
		'header' => [
			qw(Name Termini Termini Strand Upstream_Sequence SD SD_Length Spacing Start_Codon)
		],
		'data' => [],
	},
);
my $seq = $bio->getSeqIO($options->file);
my $seq_obj = $seq->next_seq;

# Use the naïve prediction algorithm that I use
use CPT::Bio::RBS;
my $rbs_predictor = CPT::Bio::RBS->new();
$rbs_predictor->set_algorithm('naive');
foreach my $feat_object ( $seq_obj->get_SeqFeatures ) {
	if (       $options->{tag} eq 'whole'
		|| $feat_object->primary_tag eq $options->{tag} )
	{
		my $name = $bio->_getIdentifier($feat_object);

		my $upstream = $bio->intelligent_get_seq($feat_object,
			upstream => $options->{lookahead_max}+10,
			parent => $seq_obj,
		);
		# Strip off actual sequence. Extraneous lc with lib update
		my $upstream_for_analysis = lc(substr($upstream, 0, $options->{lookahead_max} + 10 - $options->{lookahead_min}));
		my $upstream_missing = lc(substr($upstream, $options->{lookahead_max} + 10 - $options->{lookahead_min}, $options->{lookahead_min}));

		my @head = ($name, $feat_object->start, $feat_object->end);
		my @tail;
		if ( $feat_object->strand == "1" ) {
			push(@head, '+');
			push(@tail,
				lc($seq_obj->subseq(
						$feat_object->start,
						$feat_object->start + 2
				)),
			);
		} else {
			push(@head, '-');
			push(@tail,
				rc(
					lc(
						$seq_obj->subseq(
							$feat_object
							  ->end - 2,
							$feat_object
							  ->end
						)
					)
				),
			);
		}


		my @rbs_hits = $rbs_predictor->predict($upstream_for_analysis);
		foreach my $rbs_hit(@rbs_hits){
			push( $results{'Sheet1'}{'data'}, [
					@head,
					$rbs_hit->upstream .  $upstream_missing,
					$rbs_hit->rbs_seq,
					$rbs_hit->score,
                    $rbs_hit->separation + length($upstream_missing),
					@tail,
				]);
		}
	}
}

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'report',
	GGO => $ggo,
);
$crr_output->CRR(data => \%results);


sub rc {
	my $val = shift;
	$val = reverse($val);
	$val =~ tr/actg/qzac/;
	$val =~ tr/qz/tg/;
	return $val;
}

=head1 DESCRIPTION

Locate possible RBSs upstream of CDSs in a genome.

=cut
