#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Bio::Tools::SeqStats;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file' => 'Input file',
			{
				validate => 'File/Input',
				required => 1,
				multiple => 1,
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
	],
	'outputs' => [
		[
			'results',
			'Stop Codon Statistics',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'stop_codons',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'stopCodonStats',
		'appname' => 'Stop Codon Statistics',
		'appvers' => '1.94',
		'appdesc' =>
'Extract a table of stop codon frequencies from a Genbank file',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/multi.gbk',
			},
			outputs => {
				results => ["stop_codons.Sheet1.csv", 'test-data/outputs/stop_codons.csv'],
			}
		},
	],
);

my %types = (
	'UAG' => 'Amber',
	'UAA' => 'Ochre',
	'UGA' => 'Opal',
);

use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;
my %global_results = ();
my %global_keys = ();
foreach my $file(@{$options->{file}}){
	# For all genomes in the GBK file

	my $seqio_object = $bio->getSeqIO($file);
	while(my $seqobj = $seqio_object->next_seq()){
		my %stop_codons;
		foreach my $feat($seqobj->get_SeqFeatures()){
			if($feat->primary_tag() eq 'CDS'){
				my $seq = $feat->seq->seq();
				my $codon = substr( $seq, -3);
				$codon =~ s/T/U/g;
				$stop_codons{$codon}++;
				# Store a list of all stop codons ever used to
				# create table header
				$global_keys{$codon} = 1;
			}
		}
		$global_results{$seqobj->display_id()} = \%stop_codons;
	}
}

my @header = ('Sequence');
my @data;
my @gk = sort(keys(%global_keys));
foreach my $sc(@gk){
	push(@header, "$sc ($types{$sc})");
}
foreach my $seqid(sort(keys(%global_results))){
	my @row = ($seqid);
	my %local_scs = %{$global_results{$seqid}};
	foreach my $sc(@gk){
		# If that start codon is used in that genome
		if(defined $local_scs{$sc}){
			push(@row, $local_scs{$sc});
		}else{
			push(@row, 0);
		}
	}
	push(@data, \@row);
}

my %results = (
	'Sheet1' => {
		header => \@header,
		data => \@data,
	},
);

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'results',
	GGO => $ggo,
);
$crr_output->CRR(data => \%results);

=head1 DESCRIPTION

This tool produces a table of stop codons used in the genome, and their frequencies.

=cut
