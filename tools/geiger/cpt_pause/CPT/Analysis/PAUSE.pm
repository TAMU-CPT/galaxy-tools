package CPT::Analysis::PAUSE;

# ABSTRACT: Library for use in PAUSE analysis
use strict;
use warnings;
use Moose;
use List::Util qw(sum);
use Statistics::Descriptive;

sub max ($$) { shift; $_[ $_[0] < $_[1] ] }
sub min ($$) { shift; $_[ $_[0] > $_[1] ] }

sub derivative {
	my ($self, $data_ref) = @_;
	my @data     = @{$data_ref};
	my @new_data;
	foreach ( my $i = 0 ; $i < scalar(@data) -1 ; $i++ ) {
		$new_data[$i+1] = $data[$i+1] - $data[$i];
	}
	return \@new_data;
}

sub find_peaks {
	my ($self, %data) = @_;
	use IPC::Run3;
	use File::Temp qw/tempfile/;
	# Store to CSV File
	my @starts = @{$data{data}};
	my ($fh0, $filename0) = tempfile('galaxy.pause.XXXXXXX');
	printf $fh0 ( "%s,%s\n", 'position', 'count' );
	for (my $i=0;$i<scalar(@starts);$i++){
		printf $fh0 "%d,%d\n", $i, (defined $starts[$i] ? $starts[$i] : 0);
	}
	close($fh0);

	my ($fh, $filename) = tempfile('galaxy.pause.XXXXXXX');
	my @cmd = ('Rscript', $data{location_of_rscript_file}, $filename0, $filename, $data{snr});
	my ($in, $out, $err);
	run3 \@cmd, \$in, \$out, \$err;

	# Read in R data
	my @values;
	while(<$fh>){
		chomp;
		push(@values,$_);
	}
	close($fh);

	unlink($filename0);
	unlink($filename);
	
	return @values;
}

sub smooth {
	my ($self, $data_ref) = @_;
	my @data     = @{$data_ref};
	my @new_data;
	my $length = scalar @data;
	foreach ( my $i = 0 ; $i < $length ; $i++ ) {
		my $avg =
		  sum( @data[ $i - 20 .. $i - 1, $i + 1 .. $i + 20 ] ) / 40;
		$new_data[$i] = $avg;
	}
	return \@new_data;
}

sub histogram {
	my ( $self, %data ) = @_;

	my @coverage = @{ $data{data} };
	my @return_coverage;
	for ( my $i = 0 ; $i < scalar(@coverage) ; $i++ ) {
		my $size = $coverage[$i];
		unless ($size) { $size = 0 }
		$return_coverage[$i] = [ $i, $size, "*" x $size ];
	}
	my %results = (
		'Sheet1' => {
			headers => [qw(Base Count Plot)],
			data    => \@return_coverage,
		}
	);
	return %results;
}

sub getCoverageDensity {
	my ( $self, %data ) = @_;

	# Load the sam file
	my $sam = Bio::DB::Sam->new(
		-bam       => $data{bam},
		-fasta     => $data{genome},
		-autoindex => 1,
	);

	# Get all alignments to our indicated FASTA file
	my @alignments = $sam->get_features_by_location(
		-seq_id => $data{fasta_id},
		-start  => 1,
		-end    => $data{fasta_length}
	);

	# Set up some variables
	my $coverage_density_max_value = 0;
	my ( @coverage_density, @read_starts, @read_ends );

	# including some for statistics
	my $stat_start = Statistics::Descriptive::Sparse->new();
	my $stat_end   = Statistics::Descriptive::Sparse->new();

	# Looping over alignments
	for my $a (@alignments) {
		my $start = $a->start;
		my $end   = $a->end;

		# Increment the number of reads starting there
		$read_starts[$start]++;
		$read_ends[$end]++;

		# And increment the coverage density
		foreach ( $start .. $end ) {
			$coverage_density[$_]++;
			if ( $coverage_density[$_] >
				$coverage_density_max_value )
			{
				$coverage_density_max_value =
				  $coverage_density[$_];
			}
		}
	}
	my @start_data_for_stats;
	my @end_data_for_stats;
	for ( my $i = 0 ; $i < $data{fasta_length} ; $i++ ) {
		if ( $read_starts[$i] ) {
			push( @start_data_for_stats, $read_starts[$i] );
		}
		if ( $read_ends[$i] ) {
			push( @end_data_for_stats, $read_ends[$i] );
		}
	}
	$stat_start->add_data(@start_data_for_stats);
	$stat_end->add_data(@end_data_for_stats);

	# Lots of data to return
	use CPT::Analysis::PAUSE::ParsedSam;
	my $psam = CPT::Analysis::PAUSE::ParsedSam->new(
		coverage_density => \@coverage_density,
		read_starts      => \@read_starts,
		read_ends        => \@read_ends,
		max              => $coverage_density_max_value,
		stats_start_max  => $stat_start->max(),
		stats_end_max    => $stat_end->max(),
		stats_start_mean => $stat_start->mean(),
		stats_end_mean   => $stat_end->mean(),
		stats_start_standard_deviation =>
		  $stat_start->standard_deviation(),
		stats_end_standard_deviation => $stat_end->standard_deviation(),
	);
	return $psam;
}

no Moose;
1;
