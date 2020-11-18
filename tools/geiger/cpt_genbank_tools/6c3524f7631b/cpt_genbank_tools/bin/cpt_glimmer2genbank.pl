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
use Bio::SeqFeature::Generic;

my $libCPT  = CPT->new();
my $options = $libCPT->getOptions(
	'options' => [
		[
			'file|f',
			'Input file (fasta)',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['fasta'],
			}
		],
		[
			"glimmerfile" => "File to obtain glimmer data from",
			{
				validate => 'File/Input',
				required => 1,
			}
		],
	],
	'outputs' => [
		[
			'genbank',
			'Glimmer converted to Genbank',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'glimmer2genbank',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank',
			}
		],
	],
	'defaults' => [
		'appid'   => 'glimmer2genbank',
		'appname' => 'Glimmer To Genbank',
		'appvers' => '1.94',
		'appdesc' =>
'Merge Glimmer/Genemark with Fasta sequence to create a Genbank file',
	],
	'tests' => [
	],
);

# Load the fasta
use Bio::SeqIO;
open( my $fasta, '<', $options->{file_fasta} );
my $fasta_seqio = Bio::SeqIO->new( -fh => $fasta, -format => 'fasta' );
my ( $seq_id, $seq );
while ( my $seqobj = $fasta_seqio->next_seq() ) {
	$seq_id = $seqobj->display_id();
	$seq    = $seqobj->seq();
}
close($fasta);

# New seqio;
use Bio::SeqFeature::Generic;
my $new_seq = Bio::Seq->new( -seq => $seq, -display_id => $seq_id );

my @glimmer;
open( FH, '<', $options->{'glimmerfile'} );
while (<FH>) {
	push( @glimmer, $_ );
}
my $regex_glimmer = '^(\S+)\s+(\d+)\s+(\d+)\s+([+-])\d\s+\d+\.\d+\s*$';
foreach my $string (@glimmer) {
	my @line = split( $regex_glimmer, $string );
	my $complement = ( $line[4] eq '-' ) ? 1 : 0;
	addFeature( $new_seq, $line[2], $line[3], $line[1], $complement );

	my $feat = new Bio::SeqFeature::Generic(
		-start       => $line[2],
		-end         => $line[3],
		-strand      => ( $complement ? -1 : 1 ),
		-primary_tag => 'CDS'
	);
	$new_seq->add_SeqFeature($feat);
}

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'genbank',
	libCPT => $libCPT,
);
$crr_output->CRR(data => $new_seq);

=head1 DESCRIPTION

This tool with the output of the glimmer gene caller and the genome's fasta file, will produce a GenBank formatted file for analysis. The genbank file lacks any feature IDs, so you might want to use the CPT's gene renumbering tool to add those.

=cut
