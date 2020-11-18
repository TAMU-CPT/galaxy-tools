#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'fasta|f',
			'Input fasta file',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['fasta'],
			}
		],
		[
			'blastclust_output',
			'Output of Blastclust',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['data'],
			}
		],
	],
	'outputs' => [
		[
			'blastclust',
			'Blastclust Output',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'blastclust',
				data_format    => 'genomic/raw',
				default_format => 'fasta'
			}
		],
	],
	'defaults' => [
		'appid'   => 'BlastclustFastaCreator',
		'appname' => 'Blastclust fasta creator',
		'appvers' => '1.94',
		'appdesc' =>
		  'Create fasta files based on the groupings in blastclust',
	]
);

# Load the fasta

use Bio::SeqIO;
use CPT::Bio;
my $bio = CPT::Bio->new();
open( my $fasta, '<', $options->{fasta} );
my $fasta_seqio = $bio->getSeqIO($fasta);
my %fasta_seqs;
while ( my $seqobj = $fasta_seqio->next_seq() ) {
	$fasta_seqs{ $seqobj->display_id() } = $seqobj->seq();
}
close($fasta);

# FASTA output

my $fasta_out;

# Read in blastclust

open( my $blastclust, '<', $options->{blastclust_output} );
my $cluster_id = 0;
while (<$blastclust>) {
	chomp $_;
	my @tmp = split( /\s+/, $_ );
	foreach my $id (@tmp) {
		$fasta_out .= sprintf( ">cluster_%03d_%s\n%s\n",
			$cluster_id, $id, $fasta_seqs{$id}, );
	}
	$cluster_id++;
}
close($blastclust);

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
        name => 'blastclust',
        GGO => $ggo,
);
$output->CRR(data => $fasta_out);

=head1 DESCRIPTION

Given the fasta file input to blastclust, as well as the output of blastclust, this tool will generate a new fasta file with the sequences from the original file grouped by their name. 

For example, if you had two sequences in a cluster, C<gp01> and C<gp02>, they would be formatted as:

    >cluser_001_gp01
    ACTGACTGATCGATGCATGCAGTC
    >cluser_001_gp02
    GATGCATGCAGTCACTGACTGATC

This will allow you to use the other tools in your toolbox (e.g., Fasta-to-Tabular, "filter on lines", Tabular-to-Fasta) to select a subset of these fasta sequences to produce a fasta file with only a single cluster.

=cut
