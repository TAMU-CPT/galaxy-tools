#!/usr/bin/env perl
use strict;
use warnings;
use Carp;
use CPT::GalaxyGetOpt;

# PODNAME: phantasm_data_comparison.pl

=head1 DESCRIPTION

Tool to compare two sets of data, where each set of data has a (unique) ID column, and a data column that is of interest for analysis. For instance, given two data tables

    #ID,Val (Table N)
    A,1
    B,-4
    C,3
    D,10

    #ID,Val (Table M)
    A,-1
    B,4
    C,-3
    D,-10

This tool could accept those files for C<file_a> and C<file_b>, with C<id_col_N> = 1 and C<data_col_N> = 2. The tool would then compare these datasets against one another and produce a larger table with an NxM comparison of the two datasets.

If the user were to select the "multiplication" option for comparison, a table like the following would be generated

    #ID_N, ID_M, Val
    A,A,-1
    A,B,4
    A,C,-3
    A,D,-10
    B,A,4
    B,B,-16
    B,C,12
    B,D,40
    ...

Every comparison is included, and these results can be used in downstream analysis either in R or in the PHAnTASM squite as a comparison table.

=head2 Algorithm Choices

Some of the algorithms are possibly unfamiliar:

=over 4

=item pdiff

Percent difference is calculated as (abs(a-b)/(a+b))

=item bit_diff

In order to protect the code, all values are converted to integers and then taken an absolute value of. (-1.5 => -1 => 1).

However if you're using this option, it's generally expected you're using it for the correct reasons and know why you would use this function.

The quick test for "should you use this" is: "If I multiply every value by 3 in your dataset, will it completely screw up the data? Or will the data be fine?" If the data will be absolutely fine (with results just scaled), then you should not be using this function to analyse your data.


E.g., f(5,4) = 1 while f(5,1) = 3

=cut

my $ggo  = CPT::GalaxyGetOpt->new();
my %comps = (
	'pdiff'   => 'Calculates % difference between the two values being compared',
	'bit_diff'    => 'Calculates number of different bits.',
	'add' => 'Add the two values',
	'sub' => 'Subtract the two values (a-b)',
	'mult' => 'Multiple the two values',
	'div' => 'Divide the two values (a/b)',
	'dist' => 'Absolute distance between two numbers (abs(a-b))',
    'equal' => 'Text equality',
);
my $options = $ggo->getOptions(
	'options' => [
		[ 'file_a', 'Input Tabular Data "From"', { validate => 'File/Input', required => 1 } ],
		[ 'file_b', 'Input Tabular Data "To"', { validate => 'File/Input', required => 1 } ],
		[],
		['Data Locations'],
		[ 'id_col_a', 'Column containing row identifier in file_a. Columns numbered starting at 1', { validate => 'Int', min => '1', default => 1}],
		[ 'id_col_b', 'Column containing row identifier in file_b. Columns numbered starting at 1', { validate => 'Int', min => '1', default => 1}],
		[ 'data_col_a', 'Column containing data of interest in file_a. Columns numbered starting at 1', { validate => 'Int', min => '1', default => 2}],
		[ 'data_col_b', 'Column containing data of interest in file_b. Columns numbered starting at 1', { validate => 'Int', min => '1', default => 2}],
		[],
		['Comparison'],
		[
			"comparison_method" => "Function or method by which two data points are compared",
			{
				validate => 'Option',
				options  => \%comps,
				required => 1,
			}
		],
		['undef_val', 'Undefined value. For operations involving division, what should undefined results be set to? (e.g. 3/0 = ?). Depending on your dataset, it may be appropriate to rescale everything 1 higher (or a similar operation) to avoid zero values.', { validate => 'Float', required => 1, default => '0' }],
	],
	'outputs' => [
		[
			'comp',
			'Compared Data',
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
		'appid'   => 'ColumnarDataComparison',
		'appname' => 'Columnar Data Comparison',
		'appdesc' => 'generates a comparison map of two datasets',
		'appvers' => '1.94.0',
	],
	'tests' => [
	],
);

my ($header_a_ref, $data_a_ref) = load_data($options->{file_a}, $options->{id_col_a}, $options->{data_col_a});
my ($header_b_ref, $data_b_ref) = load_data($options->{file_b}, $options->{id_col_b}, $options->{data_col_b});

my %comp_map_subs = (
	'pdiff'    => sub { return ($_[0]+$_[1]) != 0 ? abs($_[0]-$_[1])/($_[0]+$_[1]) : 'undef' },
	'bit_diff' => sub {
				my ($q,$w) = @_;
				my $x = abs(int($q)) ^ abs(int($w));
				my $c = 0;
				if(! defined $x || $x == 0){
					return length($q) | length($w);
				}
				for(my $i = 1; $i < ($x << 1); $i <<= 1){
					if(($i & $x) > 0){
						$c++
					}
				}
				return $c;
			},
	'add'      => sub { return $_[0] + $_[1] },
	'sub'      => sub { return $_[0] - $_[1] },
	'mult'     => sub { return $_[0] * $_[1] },
	'div'      => sub { return $_[1] != 0 ? $_[0] / $_[1] : 'undef' },
	'dist'     => sub { return abs($_[0] - $_[1]) },
    'equal'    => sub { return $_[0] eq $_[1] },
);

my @result = compare(
	data_a => $data_a_ref,
	data_b => $data_b_ref,
	method => $options->{comparison_method},
);


my @header_a = @{$header_a_ref};
my @header_b = @{$header_b_ref};
my %data = (
        'Sheet1' => {
		header => [$header_a[0], $header_b[0], 'Result'],
		data => \@result,
        }
);
use CPT::OutputFiles;
my $csv_output = CPT::OutputFiles->new(
        name => 'comp',
        GGO => $ggo,
);
$csv_output->CRR(data => \%data);

sub compare {
	my (%params) = @_;
	my @results;
	foreach my $row_a(@{$params{data_a}}){
		my ($id_a, $data_a) = @{$row_a};
		foreach my $row_b(@{$params{data_b}}){
			my ($id_b, $data_b) = @{$row_b};
			push(@results,
				[
					$id_a,
					$id_b,
					$comp_map_subs{$params{method}}->($data_a, $data_b)
				],
			);

		}
	}
	return @results;
}

sub load_data{
	my ($file, $id_col, $data_col) = @_;
	open(my $fh,'<', $file);
	my (@rows, @header);
	my $row_num = 0;
	while(<$fh>){
		chomp;
		# Remove quotes
		s/"|'//g;
		# Split on tabs or commas
		my @tmp = split /\t|,/;
		if($row_num == 0){
			# Remove apostrophes
			s/"//g;
			@header = ($tmp[$id_col-1], $tmp[$data_col-1]);
			$row_num++;
			next;
		}
		# Split on tabs or commas
		push(@rows, [$tmp[$id_col-1], $tmp[$data_col-1]]);
	}
	close($fh);
	return (\@header, \@rows);
}
