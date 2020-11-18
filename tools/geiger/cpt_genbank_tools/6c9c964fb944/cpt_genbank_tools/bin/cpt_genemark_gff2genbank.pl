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
			'file_fasta' => 'Input Fasta file',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['fasta'],
			}
		],
		[
			'file_gff' => 'Input GeneMark GFF2 formatted calls',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => [ 'gff', 'gff2', 'gff3' ],
			}
		],
	],
	'outputs' => [
		[
			'genbank',
			'Merged Genbank File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'merged',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank'
			}
		],
	],
	'defaults' => [
		'appid'   => 'GenemarkFastaMergeGFF',
		'appname' => 'GenBank from GeneMarkS GFF output',
		'appvers' => '1.94',
		'appdesc' =>
		  'creates a valid genbank file from GeneMarkS and Fasta',
	],
	'tests' => [
	],
);

# Load the fasta
use Bio::SeqIO;
use CPT::Bio;
my $bio = CPT::Bio->new();
my $fasta_seqio = $bio->getSeqIO($options->{file_fasta});
my ( $seq_id, $seq );
while ( my $seqobj = $fasta_seqio->next_seq() ) {
	$seq_id = $seqobj->display_id();
	$seq    = $seqobj->seq();
}

# New seqio;
use Bio::SeqFeature::Generic;
my $new_seq = Bio::Seq->new( -seq => $seq, -display_id => $seq_id );

# Load/open GFF crap
use Bio::Tools::GFF;
open( my $fh, '<', $options->{file_gff} );
my $gffio = Bio::Tools::GFF->new( -fh => $fh, -gff_version => 2 );
my $feature;

# Loop over features
while ( $feature = $gffio->next_feature() ) {
	my $feat = new Bio::SeqFeature::Generic(
		-start       => $feature->start(),
		-end         => $feature->end(),
		-strand      => $feature->strand(),
		-primary_tag => 'CDS'
	);
	$new_seq->add_SeqFeature($feat);
}
$gffio->close();

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
        name => 'genemark',
        libCPT => $libCPT,
);
$output->CRR(data => $new_seq);

=head1 DESCRIPTION

This tool depends on the output of GeneMarkS in GFF format. With the GFF formatted GeneMarkS output and a fasta sequence, it will produce a genbank file ready for annotation. The genbank file lacks any feature IDs, so you might want to use the CPT's gene renumbering tool to add those.

=cut
