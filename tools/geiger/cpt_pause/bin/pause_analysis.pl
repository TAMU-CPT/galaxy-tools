#!/usr/bin/env perl
# PODNAME: pause_analysis.pl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;
use Bio::DB::Sam;
use File::Temp qw/ tempfile /;
use IPC::Run3 qw(run3);

my $ggo = CPT::GalaxyGetOpt->new();

=head1 DESCRIPTION

PAUSE:Analysis is a tool to run peak finding on BAM alignments. Uses genome for sequence information surrounding peaks.

This tool takes a mapped alignment of sequencing data and uses that data for CWT (Continuous Wavelet Transform) analysis in R in order to identify peaks.

=head1 PARAMETERS

=over 4

=item C<bam_sense> bam file of reads mapping to sense strand

=item C<bai_sense> bam index for C<bam_sense>

=item C<bam_antisense> bam file of reads mapping to antisense strand

=item C<bai_antisense> bam index for C<bam_antisense>

=item C<genome> fasta sequence of the genome

=item C<starts_threshold> this parameter is the signal to noise ratio for the peak detection algorithm. 10 seems to be a good value for phages, but try adjusting downwards if you are finding too few peaks. Once you've settled on a good value and found reasonable looking peaks, confer with the PAUSE:Base plot to identify your end type.

=back

=head1 KNOWN BUGS

=over 4

=item Currently the position used in the header label is incorrect. Use the number as seen in the histogram for accurate start/end position information.

=item Information is incorrect around start/end of the genome, as reads that are off either end do not map correctly. We will probably eventually add "circular genome" assumptions which will duplicate a portion of the genome on each for better mapping results.

=back

=cut

my $options = $ggo->getOptions(
	'options' => [
		[
			'bam_sense',
			'BAM aligned reads on Sense Strand',
			{ validate => 'File/Input', required => 1 }
		],
		[
			'bam_antisense',
			'BAM aligned reads on Antisense Strand',
			{ validate => 'File/Input', required => 1 }
		],
		[
			'bai_sense',
			'BAI - bam index for Sense Strand',
			{ validate => 'File/Input' }
		],
		[
			'bai_antisense',
			'BAI - bam index for Antisense Strand',
			{ validate => 'File/Input' }
		],
		[
			'genome', 'FASTA Genome',
			{ validate => 'File/Input', required => 1 }
		],
		[
			'starts_threshold',
			'SNR threshold in CWT for starts',
			{
				validate => 'Int',
				'min'    => 1,
				default  => 10
			}
		],
	],
	'outputs'  => [
		[
			'pause',
			'PAUSE Results',
			{
				validate       => 'File/Output',
				default        => 'pause',
				data_format    => 'text/html',
				default_format => 'HTML'
			}
		],
		[
			'pause_info',
			'Extra PAUSE Information',
			{
				validate       => 'File/Output',
				default        => 'pause-extra',
				data_format    => 'text/tabular',
				default_format => 'CSV'
			}
		],
	],
	'defaults' => [
		'appid'   => 'PAUSE_Base',
		'appname' => 'PAUSE:Base',
		'appdesc' =>
'Pile-up Analysis Using Starts & Ends. Reads Illumina/454 Data and attempts to determine possible genome ends based on that data.',
		'appvers' => '0.0.1',
	],
	'tests' => [
	]
);

use Bio::SeqIO;
use CPT::Analysis::PAUSE;
use CPT::Analysis::PAUSE::SVG;
my $pause = CPT::Analysis::PAUSE->new();

my $seqio_object = Bio::SeqIO->new( -file => $options->{genome}, -format => 'fasta');
while ( my $seq_object = $seqio_object->next_seq ) {
	parse_by_sequence( $seq_object->display_id(), $seq_object->length(), $seq_object->seq() );
}

sub parse_by_sequence {
	my ( $fasta_id, $fasta_length, $seq) = @_;

	my $sense_psam = $pause->getCoverageDensity(
		bam          => $options->{bam_sense},
		genome       => $options->{genome},
		fasta_id     => $fasta_id,
		fasta_length => $fasta_length,
	);

	my $antisense_psam = $pause->getCoverageDensity(
		bam          => $options->{bam_antisense},
		genome       => $options->{genome},
		fasta_id     => $fasta_id,
		fasta_length => $fasta_length,
	);
	# data of interest
	my @sense_starts       = @{ $sense_psam->read_starts() };
	my @antisense_starts   = map { $_ ? -$_ : 0 } @{ $antisense_psam->read_ends() };


	use File::Spec;
	use File::ShareDir;
	my $dir = File::ShareDir::dist_dir('CPT-Analysis-PAUSE');
	my $static_file = File::Spec->catfile($dir,'pause.Rscript');
	
	my @sense_peaks = $pause->find_peaks(
		data => \@sense_starts,
		method => 'CWT',
		snr => $options->{starts_threshold},
		location_of_rscript_file => $static_file,
	);
	my @antisense_peaks = $pause->find_peaks(
		data => \@antisense_starts,
		method => 'CWT',
		snr => $options->{starts_threshold},
		location_of_rscript_file => $static_file,
	);

	use CPT::OutputFiles;
	my $html_report = CPT::OutputFiles->new(
		name => 'pause',
		GGO => $ggo,
	);
	my $csvdata = CPT::OutputFiles->new(
		name => 'pause_info',
		GGO => $ggo,
	);

	$html_report->times_called(1);
	use CPT::Report::HTML;
	my $report = CPT::Report::HTML->new();
	$report->h1('PAUSE Analysis Results');
	foreach(@sense_peaks){
		print STDERR "$_\n";
		$report->h3(sprintf('%s (Sense)', $_));
		$report->p("Next 30 residues: " . substr($seq,$_-1,30));
		$report->a(handle($_,'sense',$seq, $sense_psam, $antisense_psam,$html_report));
	}
	foreach(@antisense_peaks){
		print STDERR "$_\n";
		$report->h3(sprintf('%s (Antisense)', $_));
		$report->p("Next 30 residues: " . reverse(substr($seq,$_-30,30)));
		$report->a(handle($_,'antisense',$seq, $sense_psam, $antisense_psam,$html_report));
	}

	# Reset filename stuff
	$html_report->times_called(0);
	$html_report->CRR(data => $report->get_content());




	# CSV Data
	my @data;
	foreach(@sense_peaks){
		push(@data, [$_, 'sense']);
	}
	foreach(@antisense_peaks){
		push(@data, [$_, 'antisense']);
	}
	my %response = (
		'Sheet1' => {
			header => [qw(Location Strand)],
			data => \@data,
		}
	);
	$csvdata->CRR(data => \%response);
}

sub handle {
	my ($idx, $strand, $seq, $sense_psam, $antisense_psam, $html_report) = @_;

	######################################33
	## Left panel contains a copy of seq + raw histo
	print STDERR "\tLeft Panel\n";
	my @starts;
	my @coverage;
	if($strand eq 'sense'){
		@starts = @{ $sense_psam->read_starts() };
		@coverage = @{ $sense_psam->coverage_density() };
	}else{
		@starts = map { $_ ? -$_ : 0 } @{ $antisense_psam->read_ends() };
		@coverage = map { -$_ } @{ $antisense_psam->coverage_density() };
	}
	my @lines;
	my $local_subseq = substr($seq,$idx-10,21);
	my $start = $idx-10;
	my $end = $start + 21;
	for(my $i = $start; $i < $end; $i++){
		my $line_add = sprintf("%-8d ", $i) . ' ';
		if($starts[$i]){
			$line_add .= '=' x abs($starts[$i]);
		}
		push(@lines,$line_add);
	}

	######################################33
	## Right panel contains a mini pause-image
	print STDERR "\tRight Panel\n";
	my $sense_start_end_max = $pause->max( $sense_psam->stats_start_max(),
		$sense_psam->stats_end_max() );
	my $antisense_start_end_max = $pause->max(
		$antisense_psam->stats_start_max(),
		$antisense_psam->stats_end_max()
	);
	my $start_end_max_num =
	  $pause->max( $sense_start_end_max, $antisense_start_end_max );
	my $max = $pause->max( $sense_psam->max(), $antisense_psam->max() );


	my $row_size           = 1000 * 10;
	my $line_height        = 200;
	my $inter_line_spacing = $line_height + 75;
	my $row_width          = 400;
	my $num_rows           = 1;
	my $x_border           = 70;
	my $y_border           = 20;

	my $psvg = CPT::Analysis::PAUSE::SVG->new(
		width  => $row_width + 2 * $x_border,
		height => $line_height * ( 1 + $num_rows ) +
		  ( $inter_line_spacing * ( $num_rows - 1 ) ) +
		  2 * $y_border,
		line_height        => $line_height,
		inter_line_spacing => $inter_line_spacing,
		vertical_offset    => 0,
		start_end_max_num  => $start_end_max_num,
		num_rows           => $num_rows,
		row_size           => $row_size,
		x_border           => $x_border,
		y_border           => $y_border,
		max                => $max,
		fasta_id           => $idx,
		row_width          => $row_width,
	);
	$psvg->setup();

	$psvg->plot_data_subset(
		from   => $start-40,
		to     => $end+40,
		regular => [
			{
				data => \@starts,
				line => 'purple',
				fill => 'none',
				name => 'Starts',
			},
		],
		rescale => [
			{
				data => \@coverage,
				line => 'black',
				fill => 'rgb(200,200,200)',
				name => 'Coverage Density',
			},
		],
	);

	$html_report->subCRR(
		filename => sprintf('%s.%s', $idx, $strand),
		extension => 'svg',
		data_format => 'image/svg',
		format_as => 'SVG',
		data => $psvg
	);


	my $svg_name = sprintf('%s.%s.svg', $idx, $strand);
	my $left_div = '<pre>' . join("\n", @lines) . '</pre>';
	my $right_div = sprintf('<object data="%s" width="540" height="540" type="image/svg+xml"><embed src="%s" width="540" height="540" type="image/svg+xml" /></object>',
		$svg_name, $svg_name);

	my $css = 'width:50%;float:left;height:544px;';
	return sprintf('<div style="%s">%s</div><div style="%s">%s</div>', $css, $left_div,$css, $right_div);
	
}

