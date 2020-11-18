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
use Data::Dumper;

my $libCPT = CPT->new();

my $options = $libCPT->getOptions(
	'options' => [
		[
			'file',
			'Input file',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
		[
			'tag',
			'Tag to extract',
			{
				required => 1,
				validate => 'Genomic/Tag',
				multiple => 1
			}
		],
		[ 'translate', 'Translate the sequence during analysis' ],
		[ 'query', 'Optional query to select features. Text entered here will be searched through all tags for selected feature type. Use ? as a single character wildcard, and * as any number of character wildcard', { validate=>'String' } ],
		[ 'n_bases_upstream', 'Adds N bases upstream to the result', { validate => 'Int', default => '0', min => '0' } ],
		[ 'n_bases_downstream', 'Adds N bases downstream to the result', { validate => 'Int', default => '0', min => '0' } ],
		[ 'strip_stops', 'Remove all stop codons from translated proteins', { validate => 'Flag'} ],

	],
	'outputs' => [
		[
			'fasta',
			'Fasta export of selected tag',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'export',
				data_format    => 'genomic/raw',
				default_format => 'Fasta',
			}
		],
	],
	'defaults' => [
		'appid'   => 'GBKFeatureExport',
		'appname' => 'Genbank Feature Export',
		'appvers' => '1.94',
		'appdesc' => 'Export features from a Genbank File',
	],
	'tests' => [
		{
			test_name    => "Default DNA export",
			params => {
				'file' => 'test-data/inputs/single.gbk',
				'tag' => 'CDS',
			},
			outputs => {
				'fasta' => ['export.fa', 'test-data/outputs/export.default.fa'],
			}
		},
		{
			test_name => "translated protein export",
			params => {
				'file' => 'test-data/inputs/single.gbk',
				'tag' => 'gene',
				'translate' => '',
			},
			outputs => {
				'fasta' => ['export.fa', 'test-data/outputs/export.translate.fa'],
			}
		},
		{
			test_name => "split features",
			params => {
				'file' => 'test-data/inputs/split-feature.gbk',
				'tag' => 'CDS',
			},
			outputs => {
				'fasta' => ['export.fa', 'test-data/outputs/export.split.fa'],
			}
		},
		{
			test_name => "split features, translated",
			params => {
				'file' => 'test-data/inputs/split-feature.gbk',
				'tag' => 'CDS',
				'translate' => '',
			},
			outputs => {
				'fasta' => ['export.fa', 'test-data/outputs/export.splittrans.fa'],
			}
		},
		{
			test_name => "Upstream/Downstream",
			params => {
				'file' => 'test-data/inputs/split-feature.gbk',
				'tag' => 'CDS',
				'translate' => '',
				'n_bases_upstream' => 10,
				'n_bases_downstream' => 50,
			},
			outputs => {
				'fasta' => ['export.fa', 'test-data/outputs/export.up_down.fa'],
			}
		},
	],
);
my ($query, $regex);
if(defined $options->{query}){
$query = $options->{query};
$query =~ s/[^A-Za-z0-9_-]*//g;
$query =~ s/\?/./g;
$query =~ s/\*/.*/g;
$regex = qr/$query/;
}

use CPT::Bio;
my $bio = CPT::Bio->new();

use Bio::SeqIO;
my $seqio_object = $bio->getSeqIO($options->{file});
# For all genomes in the GBK file
my $seqs;
while(my $seq_object = $seqio_object->next_seq){
	foreach my $feat ( $seq_object->get_SeqFeatures ) {
		if(validateFeat($feat)){
			my $header = $bio->_getIdentifier($feat);
			my $seq = $bio->intelligent_get_seq($feat,
				upstream => $options->{n_bases_upstream},
				downstream => $options->{n_bases_downstream},
				parent => $seq_object,
			);
			# Proteins come with translated stop codon
			if($options->{translate}){
				$seq = $bio->translate($seq);
				if($options->{strip_stops}){
					$seq =~ s/\*//g;
					$seq =~ s/\+//g;
					$seq =~ s/#//g;
				}
			}
			$seqs .= ">$header\n$seq\n";
		}
	}
}
use CPT::OutputFiles;
my $psmout = CPT::OutputFiles->new(
	name => 'fasta',
	libCPT => $libCPT,
);
$psmout->CRR(data => $seqs);

sub validateFeat {
	my ($feat) = @_;
	my $correct_type = $feat->primary_tag eq 'CDS';
	if(defined $options->{query}){
		foreach my $tag($feat->get_all_tags()){
			foreach my $val ($feat->get_tag_values($tag)){
				if($val =~ $regex){
					return $correct_type;
				}
			}
		}
	}else{
		return $correct_type;
	}
}

=head1 DESCRIPTION

This tool exports features from a GenBank formatted file. It has the ability to prepend and append sequence up and downstream of the feature, as well as a built in sequence translator.

Some downstream tools do not want stop codons included, and will reject sequences exported from this tool (e.g., blast). If you find this is the case, simply select C<strip_stops>

=cut
