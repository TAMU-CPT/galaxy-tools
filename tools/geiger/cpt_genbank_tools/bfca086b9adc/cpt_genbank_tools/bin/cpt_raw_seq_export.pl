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
			'file|f',
			'Input file',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
	],
	'outputs' => [
		[
			'fasta',
			'Raw Sequence',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'seq',
				data_format    => 'genomic/raw',
				default_format => 'Fasta',
			}
		],
	],
	'defaults' => [
		'appname' => 'Raw Seq Export',
		'appdi'   => 'RawSeq',
		'appvers' => '1.94',
		'appdesc' =>
		  'Export the whole genomic sequence from a Genbank file',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				"file" => "test-data/inputs/multi.gbk",
			},
			outputs => {
				"results" => ["seq.fa", "test-data/outputs/raw_export.fa"]
			},
		},
	]
);


use Bio::SeqIO;
use Bio::Seq;

my @fasta;
use CPT::Bio;
my $bio = CPT::Bio->new();
my $seqio_object = $bio->getSeqIO($options->{file});
while(my $seq_obj = $seqio_object->next_seq()){
	my $seq = $seq_obj->seq();
	my $bs = Bio::Seq->new( -seq => $seq, -display_id => $seq_obj->display_id());
	push(@fasta,$bs);
}

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'fasta',
	libCPT => $libCPT,
);
$crr_output->CRR(data => \@fasta);

=head1 DESCRIPTION

Export the entire genomic sequence from a GenBank file

=cut
