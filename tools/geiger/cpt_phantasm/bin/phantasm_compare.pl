#!/usr/bin/env perl
use strict;
use warnings;
use Carp;
use CPT::GalaxyGetOpt;

# PODNAME: phantasm_compare.pl
#
my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[ 'file', 'Input Two-Way Comparison. Should be a three-column file with column 1 and 2 being a from/to name of a node. Column 3 should be the edge weight between the nodes. For bi-directional nodes, this tool will average the values',
			{
				validate => 'File/Input',
				multiple => 1,
				required => 1,
			}
		],
		[ 'metric_weight', 'How important/how much of the overall result should this particular variable count as. If you have 3 files and you specify the numbers 2, 1, 1, then dataset 1 will have a 50% weight and the other two will have a 25% weighting.',
			{
				validate => 'Float',
				multiple => 1,
				required => 1,
				min => 1,
				default => [1],
			}
		],
		[ 'metric_multiplier', 'Dataset multiplier. Datasets are allowed to be transformed again before use with a formula of f(x) = a*x + b; this represents the value of "a"',
			{
				validate => 'Float',
				multiple => 1,
				required => 1,
				default => [1],
			}
		],
		[ 'metric_adjustment', 'Dataset constant. Datasets are allowed to be transformed again before use with a formula of f(x) = a*x + b; this represents the value of "b"',
			{
				validate => 'Float',
				multiple => 1,
				required => 1,
				default => [0],
			}
		],

	],
	'outputs' => [
		[
			'comp_map',
			'Comparison Matrix',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'comp_map',
				data_format    => 'text/tabular',
				default_format => 'TSV_U',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PHAnTASM.ComparisonMapper',
		'appname' => 'Comparison Mapper',
		'appdesc' => 'utilising multiple sets of input data as a basis for generating a comparison matrix which can then be fed into the clustering tool',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);

my @files = @{$options->{file}};
my @metric_weight = @{$options->{metric_weight}};
my @metric_multiplier = @{$options->{metric_multiplier}};
my @metric_adjustment = @{$options->{metric_adjustment}};

if (scalar @files != scalar @metric_weight || scalar @metric_weight != scalar @metric_adjustment || scalar @metric_adjustment != scalar @metric_multiplier){
	carp 'List lengths must be identical. If you want to use the default value, supply it!';
}

use List::Util qw/sum/;
use List::MoreUtils qw/each_array/;

my %result_map;

my $weight_sums = sum(@metric_weight);

# Iterate across all arrays
my $it = each_array(@files, @metric_weight, @metric_multiplier, @metric_adjustment);

while( my ($file, $weight, $mult, $adj) = $it->() ){
	my %data = load_data($file);
	foreach my $i(keys(%data)){
		foreach my $j(keys(%{$data{$i}})){
			$result_map{$i}{$j} += ($mult * $data{$i}{$j} + $adj) * ($weight/$weight_sums);
		}
	}
}

sub load_data{
	my ($file) = @_;
	open(my $fh, '<', $file);
	my $row = 0;
	my %data;
	while(<$fh>){
		chomp;
		# Ignore first line
		if($row == 0){
			$row++;
			next;
		}
		# Remove quotes if there are any
		s/"//g;
		# Split by commas or tabs
		my ($from, $to, $value) = split(/,|\t/,$_);
		$data{$from}{$to} = $value;
	}
	return %data;
}

my @keys = keys(%result_map);
my @table_main;
# Header row
my @header = ( '', @keys);
foreach my $from(@keys){
	my @row_main = ($from);
	foreach my $to(@keys){

		# Average from->to and to->from
		my $val = 0;
		my $c = 0;
		if(defined($result_map{$from}{$to})){
			$val += $result_map{$from}{$to};
			$c++;
		}
		if(defined($result_map{$to}{$from})){
			$val += $result_map{$to}{$from};
			$c++;
		}

		if($val == 0){
			$val = '0.0';
		}

		# We shouldn't encounter this. Maybe we should issue a warning?
		if($c == 0){
			$val = 'NA';
			warn "Could not find any values for $from -> $to. This is very odd. Please ensure that all of your comparisons were generated from the same source data.\n";
		}

		push(@row_main, $val);
	}
	push(@table_main, \@row_main);
}

my %data = (
        'Sheet1' => {
		header => \@header,
		data => \@table_main,
        }
);
use CPT::OutputFiles;
my $csv_output = CPT::OutputFiles->new(
	name => 'comp_map',
	GGO => $ggo,
);
$csv_output->CRR(data => \%data);

=head1 NAME

PHAnTASM Cluster Comparison Map Table Generator

=head1 DESCRIPTION

Given a set C<G> of genomes, and one or more sets of scores C<S_{1..n}>, this tool aggregates those sets of scores (with appropriate modifications) and produces a single comparison map which can be fed in to the PHAnTASM Hierarchical Clustering Tool

To put it mathematically,

    let f(S_i, g1, g2) be a function of S_i which returns the score associated with the pair g1,g2 in S_i

    forall sets of scores S_i in S_{1..n}:
        forall genomes g1,g2 in G:
            f(S_i, g1, g2) is defined

    From this knowledge, and a set of transformations T with a unique mapping to S, we can generate a result mapping R (like S)

    f(R, g1, g2) = \sum_{i=1}^{n} f(S_i, g1, g2) * T_i

    Which represents an aggregation of all our score sets S

=cut
