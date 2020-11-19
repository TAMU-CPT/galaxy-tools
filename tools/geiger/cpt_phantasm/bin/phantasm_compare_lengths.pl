#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Carp;
use CPT;
# PODNAME: phantasm_compare_lengths.pl

my $libCPT  = CPT->new();
my $options = $libCPT->getOptions(
	'options' => [
		[ 'file', 'Genome Length list',
			{
				validate => 'File/Input',
				required => 1,
			}
		],
	],
	'outputs' => [
		[
			'comp_lens',
			'Scored len pairs',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'comp_lens',
				data_format    => 'text/tabular',
				default_format => 'TSV_U',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PHAnTASM.Comparison.Lengths',
		'appname' => 'Length Scoring',
		'appdesc' => 'produces from-to-score table of lengths for input length table (CPT LOGS)',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);

my %map;
my %comp_len_map;

# Load in data
open(my $lens_fh, '<', $options->{file});
my $row = 0;
while(<$lens_fh>){
	chomp;
	# Ignore first line
	if ($row == 0){
		$row++;
		next;
	}
	# Remove quotes if there are any
	s/"//g;
	# Split by commas and tabs
	my ($seq, $len) = split(/,|\t/,$_);
	$map{$seq} = $len;
}
close($lens_fh);

my @table_data;
# Compare N vs N and save in data table
my @ks = sort(keys(%map));
foreach my $from(@ks){
	foreach my $to(@ks){
		my $score = abs($map{$from}-$map{$to})/($map{$from}+$map{$to});
		# Trick it to print a value equal to zero.
		if($score == 0){$score = '0.0';}
		# Push onto the table
		push(@table_data, [$from,$to,$score]);
	}
}
# Save Data
my %data = (
        'Sheet1' => {
		header => ['From','To', 'Score'],
		data => \@table_data,
        }
);
use CPT::OutputFiles;
my $csv_output = CPT::OutputFiles->new(
        name => 'comp_lens',
        libCPT => $libCPT,
);
$csv_output->CRR(data => \%data);

=head1 NAME

PHAnTASM Length Scoring Tool

=head1 DESCRIPTION

Given two genome lengths A and B, this tool calculates C<abs(A-B)/A+B>

=cut
