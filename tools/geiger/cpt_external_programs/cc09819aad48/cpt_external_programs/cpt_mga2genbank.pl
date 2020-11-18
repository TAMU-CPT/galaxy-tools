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
			'file_mga' => 'Input MGA formatted calls',
			{
				required    => 1,
				validate    => 'File/Input',
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
		'appid'   => 'MGAFastaMergeGenbank',
		'appname' => 'GenBank from MGA output',
		'appvers' => '1.94',
		'appdesc' =>
		  'creates a valid genbank file from MGA gene calls and Fasta',
	]
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


# Load/open GFF crap
open( my $fh, '<', $options->{file_mga} );
# Loop over features
while ( <$fh> ){
	next if ($_ =~ /^#/);
	chomp $_;
	my ($gene_id, $start, $end, $strand, $phase, $complete, $score, $model, $rbs_start, $rbs_end, $rbs_score) = split(/\t/, $_);

	if($rbs_start ne '-'){
		my $rbs = new Bio::SeqFeature::Generic(
			-start       => $rbs_start,
			-end         => $rbs_end,
			-strand      => $strand eq '+' ? 1 : -1,
			-primary_tag => 'RBS',
			-tag         => {
				locus_tag => [$gene_id],
			}
		);
		$new_seq->add_SeqFeature($rbs);
	}
	my $feat = new Bio::SeqFeature::Generic(
		-start       => $start + $phase,
		-end         => $end,
		-strand      => $strand eq '+' ? 1 : -1,
		-primary_tag => 'CDS',
		-tag         => {
			locus_tag => [$gene_id],
		}
	);
	$new_seq->add_SeqFeature($feat);
}

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
        name => 'genbank',
        libCPT => $libCPT,
);
$output->CRR(data => $new_seq);
