#!/usr/bin/perl
#
#       Code written by Eric Rasche
#	       mailto:rasche.eric@yandex.ru
#	       tel:404.692.2048
#	       http://eric.rasche.co.uk
#       for
#	       Center for Phage Technology
#

use strict;
use warnings;

use CPT;
use Bio::Seq;
use Bio::SeqIO;

my $libCPT = CPT->new();

my $options = $libCPT->getOptions(
	'options' => [
		[
			'file|f',
			'Input file',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
		[
			"tag",
			"Analyse by a specific tag",
			{ validate => 'Genomic/Tag', default => 'CDS' }
		],
		[
			'lookahead',
			'How far should the script look ahead',
			{
				required => 1,
				validate => 'Int',
				default  => 30
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
		'appvers' => '1.94',
		'appdesc' =>
'Create a table of shine-delgarno sequence upstream of proteins',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/single.gbk',
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
			qw(Name Start End Strand Upstream_Sequence SD Spacing Start_Codon)
		],
		'data' => [],
	},
);
$options->{'lookahead'}++
  ;    #Still need to add one, moved here so we have prettier names.
my $seq = Bio::SeqIO->new(
	-file   => $options->file,
	-format => 'genbank'
);
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
			upstream => $options->{lookahead},
			parent => $seq_obj,
		);
		# Strip off actual sequence. Extraneous lc with lib update
		$upstream = lc(substr($upstream, 0, $options->{lookahead}));

		my @SDs = $rbs_predictor->predict($upstream);

		my @tmp;
		if ( $feat_object->strand == "1" ) {
			@tmp = (
				$name,    #name
				$feat_object->start,
				$feat_object->end,
				"+", @SDs,
				lc(
					$seq_obj->subseq(
						$feat_object->start,
						$feat_object->start + 2
					)
				),
			);
		}
		else {
			@tmp = (
				$name,    #name
				$feat_object->end,
				$feat_object->start,
				"-", @SDs,
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

		push( $results{'Sheet1'}{'data'}, \@tmp );
	}
}

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'report',
	libCPT => $libCPT,
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
