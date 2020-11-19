#!/usr/bin/env perl
use strict;
use warnings;
use Carp;
use CPT::GalaxyGetOpt;
use Algorithm::KMeans;
use List::MoreUtils qw/minmax/;

# PODNAME: phantasm_data_rescale.pl

=head1 DESCRIPTION

PHAnTASM Data Rescaling tool will help the user fix the distribution of data in their tables of values. Given a table

    #ID, VAL
    A,-1
    B,-2
    C,3
    D,4

this tool will allow you to:

=over 4

=item Change scale of data

=item Cluster data

=back

=head2 Changing Data Scale

Rescaling data allows you to adapt your data into a format more paletable to downstream analysis. Given a set of negative numbers, it may be desirable to shift these to be positive numbers so a logarithm may be applied. By setting C<min> and C<max>, it is possible to do exactly that. For the above table, if we set C<min = 1> and C<max = 10>, we would see the following results:

    A,2.5
    B,1
    C,8.5
    D,10

Very easy!

=head2 Clustering Data

You may find that the range of values in your dataset is too larger, and that these values should be binned or clustered. This tool offers a couple different methods by which that might be achieved.

=over 4

=item User specified breakpoints

These are the simplest and most controlled option. The user specifies a number of different breakpoints and values will be distributed into bins that are inside those breakpoints.

=item % based paritions

Another very simple option is to say "I want N bins". The number line along which your data is distributed is then split into N+1 pieces and your data is binned accordingly.

=item k-means clustering

Alternatively, you can use k-means clustering to handle your analysis. K-means clustering will attempt to cluster your data into K bins, which you can specify or let the library determine automatically. If you specify K, be sure that

    N > 2 * K^2

where N is the number of data points you have and K is the number of clusters you wish to have.

=cut

my %part_type = (
	'none' => 'No partitioning is done',
	'user' => 'User specified breakpoints',
	'kmeans' => '1D K-means clustering',
	'simple' => 'Simple % of min/max based partitions',
);

my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[ 'file', 'Input Tabular Data', { validate => 'File/Input', required => 1 } ],
		[ 'column', 'Selected column from dataset to apply transformation to, where 1 is the first column', { validate => 'Int', min => 1, required => 1, default => 2}],
		[],
		['Operations'],
		['min', 'New minimum value for data set. All data will be shifted upwards by a factor calculated from the data set\'s max valu', { validate => 'Int'} ],
		['max', 'New maximum value for data set. All data will be shifted downward by a factor calculated from the data set\'s min value', { validate => 'Int'} ],
		[],
		['Partitioning Opertions'],
		['none', 'No partitioning is done'],
		['partition_type', 'Select the type of partitioning to use (galaxy only)', { validate => 'Option', options => \%part_type }],
		['user_provided_breakpoints', 'This option allows you to specify a set of values that act as breakpoints for the data. Any value less than OR equal to will be within the cutoff. Data is then renumbered as 0..#_of_clusters. If these numbers should be different, re-run this tool on the output with min/max/invert, or use the PHAnTASM Data Transform.', { validate => 'Float', multiple => 1 }],
		['k_means_clusters', 'Number of k-means clusters to create (or 0 for "best guess")', { validate => 'Int' }],
		['percentage_based_partitions', 'Simple partitions based on an N equal subdivisions of min/max. Good for evenly distributed data.' , { validate=>'Int'}],
	],
	'outputs' => [
		[
			'rescale',
			'Rescaled Data',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'rescaled',
				data_format    => 'text/tabular',
				default_format => 'TSV_U',
			}
		],
	],
	'defaults' => [
		'appid'   => 'ColumnarDataRescale',
		'appname' => 'Columnar Data Rescale and Partitioning',
		'appdesc' => 'transforms data sets according to a particular modification type to coerce them for use in downstream analysis',
		'appvers' => '1.94.2',
	],
	'tests' => [
	],
);

my ($header_ref, $data_ref) = load_data($options->{file});
my $result = partition(
	data => $data_ref,
	basic_transforms => {
		min => $options->{min},#int
		max => $options->{max},#int
		log => $options->{log}, #bool
		invert => $options->{invert}, #bool
	},
	partitions => {
		user_provided_breakpoints => $options->{user_provided_breakpoints},
		k_means_clusters => $options->{k_means_clusters},
		percentage_based_partitions => $options->{percentage_based_partitions},
	},
);


my %data = (
        'Sheet1' => {
		header => $header_ref,
                data => $result,
        }
);
use CPT::OutputFiles;
my $csv_output = CPT::OutputFiles->new(
        name => 'rescale',
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

sub partition{
	my (%params) = @_;
	#my @data = @{$params{data}};
	my %basic = %{$params{basic_transforms}};
	my %partitions = %{$params{partitions}};

	my ($squash, $offset) = calculate_reshape(
		min => $basic{min},
		max => $basic{max},
		data => $params{data},
	);
	apply_reshape(
		data => $params{data},
		squash => $squash,
		offset => $offset,
	);
	

	partition_data(
		data => $params{data},
		user => $partitions{user_provided_breakpoints},
		kmeans => $partitions{k_means_clusters},
		percent => $partitions{percentage_based_partitions},
	);
	return $params{data};
}

sub partition_data {
	my (%params) = @_;
	if(defined($params{user})){
		parition_breakpoints(
			data => $params{data},
			user => $params{user},
		);
	}
	if(defined($params{kmeans})){
		require File::Temp;
		my ($tmp_fh, $tmp_fn) = File::Temp::tempfile();
		# Write data to file
		my @revmap;
		my $idx = 0;
		foreach my $row(@{$params{data}}){
			my @rowdata = @{$row};
			printf $tmp_fh "%s\t%s\n", $idx, $rowdata[$options->{column} - 1];
			# Store an idx to row mapping so we can recover this data later and updte our array
			$revmap[ $idx++ ] = $row;
		}
		close($tmp_fh);

		# Run the clustering
		my $clusterer = Algorithm::KMeans->new(
			datafile => $tmp_fn,
			mask => "N1",
			K => $params{kmeans},
			terminal_output => 0,
		);
		$clusterer->read_data_from_file();
		$clusterer->kmeans();
		my ($clusters, $cluster_centers) = $clusterer->kmeans();
		my %cluster_center_map;
		for(my $j = 0; $j < scalar(@{$clusters}); $j++){
			$cluster_center_map{${$cluster_centers}[$j][0]} = ${$clusters}[$j];
		}

		# Apply cluster numbers
		my $cluster_num = 0;
		foreach my $cluster_center(sort{$a<=>$b} keys(%cluster_center_map)){
		#foreach my $cluster(@{$clusters}){
			# This will return an array of the row idxs we created above
			my $cluster = $cluster_center_map{$cluster_center};
			foreach my $value(@{$cluster}){
				my @row = @{$revmap[ $value ]};
				$row[$options->{column} - 1] = $cluster_num;
				$revmap[ $value ] = \@row;
			}
			$cluster_num++;
		}

		# Copy back to the data
		my $idx2 = 0;
		foreach my $row(@{$params{data}}){
			$row = $revmap[$idx2++];
		}

		return;
	}
	if(defined($params{percent})){
		my ($data_min, $data_max) = get_min_max(data => $params{data});
		# Correct for smallest item always being binned separately.
		$data_min -= 1;
		my $width = abs($data_max - $data_min) / $params{percent};
		my @breakpoints;
		for(my $i = $data_min; $i < $data_max; $i+= $width){
			push(@breakpoints, $i);
		}
		use Data::Dumper;
		print Dumper \@breakpoints;

		parition_breakpoints(
			data => $params{data},
			user => \@breakpoints,
		);

		return;
	}
	return;
}

sub parition_breakpoints {
	my (%params) = @_;

	my %breakpoint_map;
	# Set up breakpoint array
	my $last_break = '-inf';
	my $i = 0;
	foreach my $break(sort {$a <=> $b} @{$params{user}}){
		$breakpoint_map{$last_break}{$break} = $i++;
		$last_break = $break;
	}
	$breakpoint_map{$last_break}{'+inf'} = $i++;

	foreach my $row(@{$params{data}}){
		my @rowdata = @{$row};
		my $val = $rowdata[$options->{column} - 1];
		my $new_val;

		my $hit = 0;
		foreach my $left_side(keys(%breakpoint_map)){
			foreach my $right_side(keys($breakpoint_map{$left_side})){
				if($left_side eq '-inf'){
					if($val <= $right_side){
						$hit = 1;
						$new_val = $breakpoint_map{$left_side}{$right_side};
					}
				}
				elsif($right_side eq '+inf'){
					if($val > $left_side){
						$hit = 1;
						$new_val = $breakpoint_map{$left_side}{$right_side};
					}
				}elsif($val <= $right_side && $val > $left_side){
					$hit = 1;
					$new_val = $breakpoint_map{$left_side}{$right_side};
				}

				if($hit){
				}
			}
		}
		if(!$hit){
			$val = 'bad';
		}else{
			$val = $new_val;
		}
		$rowdata[$options->{column} - 1] = $val;
		$row = \@rowdata;
	}
}

sub apply_reshape {
	my (%params) = @_;
	printf "Squashing with %s and %s\n", $params{squash}, $params{offset};
	foreach my $row(@{$params{data}}){
		my @rowdata = @{$row};
		$rowdata[$options->{column} - 1] = ($params{squash} * $rowdata[$options->{column} - 1]) + $params{offset};
		$row = \@rowdata;
	}
}

sub get_min_max {
	my (%params) = @_;
	my @data;
	foreach my $row(@{$params{data}}){
		push(@data, ${$row}[$options->{column} - 1]);
	}
	return minmax(@data);
}

sub calculate_reshape {
	my (%params) = @_;
	# If we have a new min
	my ($data_min, $data_max) = get_min_max(data => $params{data});
	my $new_min = $data_min;
	if(defined $params{min}){
		$new_min = $params{min};
	}
	# If we have a new max
	my $new_max = $data_max;
	if(defined $params{max}){
		$new_max = $params{max};
	}
	# Calculate the scaling (squash)
	my $squash = abs($new_max - $new_min) / abs($data_max - $data_min);

	# Calculate the offset
	my $offset = 0;
	if(defined $params{min}){
		$offset += $params{min} - $squash * $data_min;
	}elsif(defined $params{max}){
		$offset += $params{max} - $squash * $data_max;
	}
	return ($squash, $offset);
}
