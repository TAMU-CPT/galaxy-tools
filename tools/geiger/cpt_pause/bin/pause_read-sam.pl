#!/usr/bin/env perl
# PODNAME: pause_read-sam.pl
use strict;
use Bio::DB::Sam;
use CPT::GalaxyGetOpt;
use Data::Dumper;
use Statistics::Descriptive;

my $ggo = CPT::GalaxyGetOpt->new();

=head1 DESCRIPTION

PAUSE:Plot is a tool to generate static visualizations of read alignments to a genome. No automated analysis is done with this tool.

=head1 PARAMETERS

=over 4

=item C<bam_sense> bam file of reads mapping to sense strand

=item C<bam_antisense> bam file of reads mapping to antisense strand

=item C<genome> fasta sequence of the genome

=item C<starts_threshold> this parameter is the signal to noise ratio for the peak detection algorithm. 10 seems to be a good value for phages, but try adjusting downwards if you are finding too few peaks. Once you've settled on a good value and found reasonable looking peaks, confer with the PAUSE:Base plot to identify your end type.

=back

=head2 PLOT OPTIONS

=over 4

=item C<kb_per_row> This parameter controls segmentation of the genome as plots are cut into separate rows for easier plotting/viz.

=item C<line_height> Size of each row of the genome

=item C<plot_width> Controls overall width of the plot. Can allow for smaller/bigger plots as needed.

=item C<x_border> Amount of blank space on the left and right sides of the plot.

=item C<y_border> Amount of blank space on the top and bottom sides of the plot.

=item C<nopng> Do not generate a PNG (only an SVG)

=back

=cut

my $options = $ggo->getOptions(
	'options' => [
		[
			'fasta',
			'Input genome (FASTA)',
			{ required => 1, validate => 'File/Input' }
		],
		[
			'sense_bam',
			'Sense BAM File',
			{ required => 1, validate => 'File/Input' }
		],
		[
			'antisense_bam',
			'Antisense BAM File',
			{ required => 1, validate => 'File/Input' }
		],
		[],
		["Plot Options"],
		['kb_per_row'  , 'Base number of kb plotted per row'           , { required => 1                                           , validate => 'Int' , default => 20 }]  ,
		['line_height' , 'Vertical size of each section of the genome' , { required => 1                                           , validate => 'Int' , default => 200 }] ,
		['plot_width'  , 'Width of plot                                , roughly in pixels (not including scales and axis labels)' , { required => 1   , validate => 'Int' , default => 1500}] ,
		['x_border'    , 'Spacing on left and right sides'             , { required => 1                                           , validate => 'Int' , default => 100}]  ,
		['y_border'    , 'Spacing on top and bottom'                   , { required => 1                                           , validate => 'Int' , default => 150}]  ,
		['nopng'       , 'Do not automatically convert SVG to PNG']    ,
	],
	'outputs'  => [
		[
			'plot',
			'Output Plot',
			{
				validate       => 'File/Output',
				default        => 'pause-plot',
				data_format    => 'text/html',
				default_format => 'HTML'
			}
		],
	],
	'defaults' => [
		'appid'   => 'read_plotter',
		'appname' => 'Aligned BAM PAUSE Plotter',
		'appdesc' =>
		  'plot data to allow visual verification of PAUSE results',
	],
	'tests' => [
	]
);

use CPT::Analysis::PAUSE;
my $pause = CPT::Analysis::PAUSE->new();

use CPT::Analysis::PAUSE::SVG;
use Bio::SeqIO;
my $seqio_object = Bio::SeqIO->new( -file => $options->{fasta} , -format => 'fasta');
while ( my $seq_object = $seqio_object->next_seq ) {
	parse_by_sequence( $seq_object->display_id(), $seq_object->length() );
}

sub parse_by_sequence {
	my ( $fasta_id, $fasta_length ) = @_;

	my $sense_psam = $pause->getCoverageDensity(
		bam          => $options->{sense_bam},
		genome       => $options->{genome},
		fasta_id     => $fasta_id,
		fasta_length => $fasta_length,
	);

	my $antisense_psam = $pause->getCoverageDensity(
		bam          => $options->{antisense_bam},
		genome       => $options->{genome},
		fasta_id     => $fasta_id,
		fasta_length => $fasta_length,
	);
	my $result_max =
	  $pause->max( $sense_psam->max(), $antisense_psam->max() );
	my @sense_coverage = @{ $sense_psam->coverage_density() };
	my @antisense_coverage =
	  map { -$_ } @{ $antisense_psam->coverage_density() };

	#@sense_coverage = @{smooth(\@sense_coverage)};
	#@antisense_coverage = @{smooth(\@antisense_coverage)};
	my $max = $pause->max( $sense_psam->max(), $antisense_psam->max() );

	my @sense_starts     = @{ $sense_psam->read_starts() };
	my @sense_ends       = @{ $sense_psam->read_ends() };
	# read_end/read_start are switched here because the backing library
	# doesn't know or care about something being on the opposite strand
	my @antisense_starts = map { -$_ } @{ $antisense_psam->read_ends() };
	my @antisense_ends   = map { -$_ } @{ $antisense_psam->read_starts() };


	use POSIX qw(ceil);
	my $row_size           = 1000 * $options->{kb_per_row};
	my $line_height        = $options->{line_height};
	my $inter_line_spacing = $line_height + 75;
	my $row_width          = $options->{plot_width};
	my $num_rows           = int( ceil( $fasta_length / $row_size ) );
	my $x_border           = $options->{x_border};
	my $y_border           = $options->{y_border};

	my $sense_start_end_max = $pause->max( $sense_psam->stats_start_max(),
		$sense_psam->stats_end_max() );
	my $antisense_start_end_max = $pause->max(
		$antisense_psam->stats_start_max(),
		$antisense_psam->stats_end_max()
	);

	my $start_end_max_num =
	  $pause->max( $sense_start_end_max, $antisense_start_end_max );

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
		fasta_id           => $fasta_id,
		row_width          => $row_width,
	);
	$psvg->setup();

	# 10kb per 500 pixels (our width)
	# create an SVG object

	$psvg->plot_data(
		regular => [
			{
				data => \@sense_starts,
				line => 'red',
				fill => 'none',
				name => 'Sense Starts',
			},
			{
				data => \@sense_ends,
				line => 'green',
				fill => 'none',
				name => 'Sense Ends',
			},
			{
				data => \@antisense_starts,
				line => 'rgb(255,0,255)',
				fill => 'none',
				name => 'Antisense Starts',
			},
			{
				data => \@antisense_ends,
				line => 'blue',
				fill => 'none',
				name => 'Antisense Ends',
			},
		],
		rescale => [
			{
				data => \@sense_coverage,
				line => 'black',
				fill => 'rgb(100,100,100)',
				name => 'Sense Coverage Density',
			},
			{
				data => \@antisense_coverage,
				line => 'black',
				fill => 'rgb(200,200,200)',
				name => 'Antisense Coverage Density',
			},
		],
	);
	use CPT::OutputFiles;
	my $html_output = CPT::OutputFiles->new(
		name => 'plot',
		GGO => $ggo,
	);

	my $html_page = '<html><head></head><body><a href="pause.svg">SVG Version</a> (right click + "Save Link As" to download)<img src="pause.png" /></body></html>';
	$html_output->CRR(
		data => $html_page,
	);

	# Output and retrieve the filename
	my ($file_location) = $html_output->subCRR(
		filename => 'pause',
		extension => 'svg',
		data => $psvg,
		data_format => 'image/svg',
		format_as => 'SVG',
	);


	# Convert the image to PNG to display
	if(!$options->{nopng}){
		my ($file_location_png) = $html_output->subCRR(
			filename => 'pause',
			extension => 'png',
			data => "",
			data_format => 'Dummy',
			format_as => 'Dummy',
		);
		qx{convert $file_location $file_location_png};
	}
}
