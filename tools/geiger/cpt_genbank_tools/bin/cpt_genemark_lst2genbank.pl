#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
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
			'file_lst' => 'Input genemark LST formatted calls',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
	],
	'outputs' => [
		[
			'genemark',
			'Merged Genbank File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'merged.gbk',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank'
			}
		],
	],
	'defaults' => [
		'appid'   => 'GenemarkFastaMergeLST',
		'appname' => 'GenBank from GeneMarkS LST output',
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

# Open LST File
my @genemark;
open( FH, '<', $options->{'file_lst'} );
while (<FH>) {
	push( @genemark, $_ );
}

my $regex_genemark =
qr/^\s*(?<gene>\d+)\s+(?<strand>[+-])\s+(?<leftextra>[<>]*)(?<left>\d+)\s+(?<rightextra>[<>]*)(?<right>\d+)\s+(?<len>\d+)\s+(?<class>\d)\s*$/;
foreach my $string (@genemark) {
	if ( $string =~ $regex_genemark ) {

		#my @line = split($regex_genemark,$string);
		my $feat = new Bio::SeqFeature::Generic(
			-start       => $+{left},
			-end         => $+{right},
			-strand      => ( $+{strand} eq '-' ? -1 : 1 ),
			-primary_tag => 'CDS'
		);
		$new_seq->add_SeqFeature($feat);
	}
}

# Write out our sequence
use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'genemark',
	GGO => $ggo,
);
$crr_output->CRR(data => $new_seq);
=head1 DESCRIPTION

This tool depends on the output of GeneMarkS in LST format. With the LST formatted GeneMarkS output and a fasta sequence, it will produce a genbank file ready for annotation. The genbank file lacks any feature IDs, so you might want to use the CPT's gene renumbering tool to add those.

=cut