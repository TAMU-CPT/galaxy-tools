#!/usr/bin/perl
# PODNAME: cpt_tmhmm_post_process.pl
use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'fasta' => 'Input file',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['fasta'],
			}
		],
		[
			'file' => 'TMHMM txt output file',
			{
				required => 1,
				validate => 'File/Input',
			}
		],
		[
			'min_length',
			'Minimum Length',
			{ validate => 'Int', required => 1, default => 10 }
		],
		[
			'max_length',
			'Maximum Length',
			{ validate => 'Int', required => 1, default => 200 }
		],
		[
			'min_exp_aas',
			'Minimum Exp number of AAs in TMHs',
			{ validate => 'Int', required => 1, default => 1 }
		],
		[
			'max_exp_aas',
			'Maximum Exp number of AAs in TMHs',
			{ validate => 'Int', required => 1, default => 50 }
		],
		[
			'min_total_prob_n_in',
			'Minimum Total prob of N-in',
			{ validate => 'Int', required => 1, default => 0 }
		],
		[
			'max_total_prob_n_in',
			'Maximum Total prob of N-in',
			{ validate => 'Int', required => 1, default => 100}
		],
		[
			'min_exp_num_first_60',
			'Minimum Exp number, first 60 AAs',
			{ validate => 'Int', required => 1, default => 0}
		],
		[
			'max_exp_num_first_60',
			'Maximum Exp number, first 60 AAs',
			{ validate => 'Int', required => 1, default => 100}
		],
		[
			'min_num_predict_tmhs',
			'Minimum Number of predicted TMHs',
			{ validate => 'Int', required => 1, default => 1}
		],
		[
			'max_num_predict_tmhs',
			'Maximum Number of predicted TMHs',
			{ validate => 'Int', required => 1, default => 10}
		],
	],
	'outputs' => [
		[
			'results',
			'Analysis results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'processed_tmhmm',
				data_format    => 'genomic/raw',
				default_format => 'Fasta'
			}
		],
		[
			'last_ones',
			'Last ones (?)',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'processed_tmhmm_lastones',
				data_format    => 'genomic/raw',
				default_format => 'Fasta'
			}
		],
	],
	'defaults' => [
		'appid'   => 'tmhmm_post_process',
		'appname' => 'TMHMM+Fasta processing',
		'appdesc' => 'uses TMHMM output and original fasta sequence to process out various portions of the TMHMM results for use later',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);

use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;
my $seqio = $bio->getSeqIO($options->{fasta});
my %fasta;
while ( my $seqobj = $seqio->next_seq() ) {
	$fasta{ $seqobj->display_id() } = $seqobj->seq();
}

my %data;
open( my $fh, '<', $options->{file} );
while (<$fh>) {
	chomp $_;
	# Comment but might be useful?
	if ( $_ =~ /^#/ ) {
		if ( $_ =~ /# ([^ ]+) ([^:]+):\s*([0-9.e-]+)/ ) {
			$data{$1}{extra}{$2} = $3;
		}
	}
	else {
		if ( $_ =~
/([^\t]*)\tTMHMM2\.0\t(inside|TMhelix|outside)\s+(\d+)\s+(\d+)\s*$/
		  )
		{
			push(
				@{ $data{$1}{hit} },
				{
					id    => $1,
					type  => $2,
					start => $3,
					end   => $4,
				}
			);
		}
	}
}
close($fh);

my %filtered_data = ();
foreach ( keys(%data) ) {

	#len
	if ( $options->{min_length} ) {
		if ( $data{$_}{extra}{'Length'} < $options->{min_length} ) {
			next;
		}
	}
	if ( $options->{max_length} ) {
		if ( $data{$_}{extra}{'Length'} > $options->{max_length} ) {
			next;
		}
	}

	# exp aas
	if ( $options->{min_exp_aas} ) {
		if ( $data{$_}{extra}{'Exp number of AAs in TMHs'} <
			$options->{min_exp_aas} )
		{
			next;
		}
	}
	if ( $options->{max_exp_aas} ) {
		if ( $data{$_}{extra}{'Exp number of AAs in TMHs'} >
			$options->{max_exp_aas} )
		{
			next;
		}
	}

	# tot prob
	if ( $options->{min_total_prob_n_in} ) {
		if ( $data{$_}{extra}{'Total prob of N-in'} <
			$options->{min_total_prob_n_in} )
		{
			next;
		}
	}
	if ( $options->{max_total_prob_n_in} ) {
		if ( $data{$_}{extra}{'Total prob of N-in'} >
			$options->{max_total_prob_n_in} )
		{
			next;
		}
	}

	# exp num first 60
	if ( $options->{min_exp_num_first_60} ) {
		if ( $data{$_}{extra}{'Exp number, first 60 AAs'} <
			$options->{min_exp_num_first_60} )
		{
			next;
		}
	}
	if ( $options->{max_exp_num_first_60} ) {
		if ( $data{$_}{extra}{'Exp number, first 60 AAs'} >
			$options->{max_exp_num_first_60} )
		{
			next;
		}
	}

	# num pred tmhms
	if ( $options->{min_num_predict_tmhs} ) {
		if ( $data{$_}{extra}{'Number of predicted TMHs'} <
			$options->{min_num_predict_tmhs} )
		{
			next;
		}
	}
	if ( $options->{max_num_predict_tmhs} ) {
		if ( $data{$_}{extra}{'Number of predicted TMHs'} >
			$options->{max_num_predict_tmhs} )
		{
			next;
		}
	}

	$filtered_data{$_} = $data{$_};
}

my $result = "";
my $last_ones;
foreach ( keys(%filtered_data) ) {

	# Construct final fasta sequences
	my $idx         = 0;
	my $last_result = "";

	my $gene_seq = $fasta{$_};
	for my $item ( @{ $data{$_}{hit} } ) {
		my %i = %{$item};
		$idx++;
		my $substr =
		  substr( $gene_seq, $i{start} - 1, $i{end} - $i{start} + 1 );
		my $str = sprintf(
			">%s_%s_%s [%s %s]\n%s\n",
			$_, $i{type}, $idx, $i{start} - 1,
			$i{end}, $substr,
		);
		$last_result = $str;
		$result .= $str;
	}
	$last_ones .= $last_result;
}
use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'last_ones',
	GGO => $ggo,
);
$crr_output->CRR(data => $last_ones);

$crr_output = CPT::OutputFiles->new(
	name => 'results',
	GGO => $ggo,
);
$crr_output->CRR(data => $result);
