#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'fasta' => 'Input fasta file',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['fasta'],
			}
		],
		[
			"file" => "SignalP txt output file",
			{
				validate => 'File/Input',
				required => 1,
			}
		],
	],
	'outputs' => [
		[
			'results',
			'Processed Sequences',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'processed_signalp',
				data_format    => 'genomic/raw',
				default_format => 'Fasta'
			}
		],
	],
	'defaults' => [
		'appid'   => 'signalp_post_process',
		'appname' => 'SignalP+Fasta processing',
		'appdesc' => 'uses SignalP output and original fasta sequence to process out various portions of the SignalP results for use later',
	],
	'tests' => [
	],
);

use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;
my $seqio = $bio->getSeqIO($options->{fasta});
my %fasta; while( my $seqobj = $seqio->next_seq()){
	$fasta{$seqobj->display_id()} = $seqobj->seq(); 
}

my $result = "";
open(my $fh, '<', $options->{file});
while(<$fh>){
	chomp $_;
	if($_ =~ /Name=([^ ]+)\s+SP='YES' Cleavage site between pos. (\d+) and (\d+):/){
		my $seq = $fasta{$1};
		$result .= sprintf(">%s_first\n%s\n", $1, substr($seq,0, $2-1));
		$result .= sprintf(">%s_second\n%s\n", $1, substr($seq,$2-1));
	}
}
close($fh);

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'results',
	GGO => $ggo,
);
$crr_output->CRR(data => $result);
