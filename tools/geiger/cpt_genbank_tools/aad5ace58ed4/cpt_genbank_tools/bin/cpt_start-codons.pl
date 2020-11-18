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
			'Start Codon Statistics',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'start_codons',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'StartCodonStats',
		'appname' => 'Start Codon Statistics',
		'appvers' => '1.94',
		'appdesc' =>
'Extract a table of start codon frequencies from a Genbank file',
	],
	'tests' => [
		{
			test_name    => 'Default',
			params => {
				'file' => 'test-data/inputs/multi.gbk',
			},
			outputs => {
				results => ['start_codons.Sheet1.csv', 'test-data/outputs/start_codons.csv'],
			}
		},
	],
);
use CPT::Bio;
my $bio = CPT::Bio->new();

use Bio::SeqIO;
my @codon_usage_data;
my %global_results = ();
my %global_keys = ();
foreach my $file(@{$options->{file}}){
	# For all genomes in the GBK file

	my $seqio_object = $bio->getSeqIO($file);
	while(my $seqobj = $seqio_object->next_seq()){
		my %start_codons;
		foreach my $feat($seqobj->get_SeqFeatures()){
			if($feat->primary_tag() eq 'CDS'){
				my $codon = substr( $feat->seq->seq(), 0, 3 );
				$start_codons{$codon}++;
				# Store a list of all start codons ever used to
				# create table header
				$global_keys{$codon} = 1;
			}
		}
		$global_results{$seqobj->display_id()} = \%start_codons;
	}
}

my @header = ('Sequence');
my @data;
my @gk = sort(keys(%global_keys));
foreach my $sc(@gk){
	push(@header, $sc);
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

This tool produces a table of start codons used in a genome

=cut
