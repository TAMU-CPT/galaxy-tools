#!/usr/bin/env perl
use strict;
use warnings;
use Carp;
use CPT::GalaxyGetOpt;

# PODNAME: phantasm_data_transform.pl

=head1 DESCRIPTION

PHAnTASM Data Transformation tool will apply simple mathematical transforms to your data. Additionally, this tool allows for chaining transforms. They will be applied sequentially, in the order that they are specified.

=cut

my $ggo  = CPT::GalaxyGetOpt->new();
my %transformations = (
	'none'   => 'x',
	'neg'    => '-x',
	'log'    => 'log(x)',
	'inv'    => '1/x',
	'abs'    => 'abs',
);
my $options = $ggo->getOptions(
	'options' => [
		[ 'file', 'Input Tabular Data', { validate => 'File/Input', required => 1 } ],
		[ 'column', 'Selected column from dataset to apply transformation to, where 1 is the first column', { validate => 'Int', min => 1, required => 1, default => 2}],
		[
			"data_transform" => "Transformation to apply to this dataset. Transformations will be applied in the order they're received. For instance, if you apply the '-x' transformation followed by 'log(x)' transformation, it will evaluate to 'log(-x)'",
			{
				validate => 'Option',
				options  => \%transformations,
				multiple => 1,
				required => 1,
			}
		],
	],
	'outputs' => [
		[
			'transformed',
			'Transformed Data',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'transformed_data',
				data_format    => 'text/tabular',
				default_format => 'TSV_U',
			}
		],
	],
	'defaults' => [
		'appid'   => 'ColumnarDataTransformation',
		'appname' => 'Columnar Data Transformation',
		'appdesc' => 'apply mathetmatical transformations to datasets',
		'appvers' => '1.94.2',
	],
	'tests' => [
	],
);

my ($header_ref, $data_ref) = load_data($options->{file});
my @result = transform(
	data => $data_ref,
	transformations => $options->{data_transform},
	column => $options->{column},
);


my %data = (
        'Sheet1' => {
		header => $header_ref,
                data => \@result,
        }
);
use CPT::OutputFiles;
my $csv_output = CPT::OutputFiles->new(
        name => 'transformed',
        GGO => $ggo,
);
$csv_output->CRR(data => \%data);



sub load_data{
	my ($file) = @_;
	open(my $fh,'<', $file);
	my (@rows, @header);
	my $row_num = 0;
	while(<$fh>){
		chomp;
		if($row_num == 0){
			# Remove apostrophes
			s/"//g;
			# Split on tabs or commas
			@header = split /\t|,/;
			$row_num++;
			next;
		}
		# Remove apostrophes
		s/"//g;
		# Split on tabs or commas
		my (@row) = split /\t|,/;
		push(@rows, \@row);
	}
	close($fh);
	return (\@header, \@rows);
}

sub transform {
	my (%params) = @_;
	my @data = @{$params{data}};
	my @transformations = @{$params{transformations}};
	my $column = $params{column};

	foreach my $row(@data){
		my @rowdata = @{$row};
		
		my $val = $rowdata[$column - 1]; # 1 indexed
		foreach my $transform(@transformations){
			$val = apply_transform($transform, $val);
		}

		$rowdata[$column - 1] = $val;
		$row = \@rowdata;
	}
	return @data;
}

sub apply_transform {
	my ($trans, $val) = @_;
	if($trans eq 'none'){
		return $val;
	}elsif($trans eq 'neg'){
		return -$val;
	}elsif($trans eq 'log'){
		return log($val);
	}elsif($trans eq 'inv'){
		return 1/$val;
	}else{
		carp 'Unknown transformation';
	}
}
