#!/usr/bin/perl
#
#       Code written by Eric Rasche
#               mailto:rasche.eric@yandex.ru
#               tel:   404.692.2048
#               http://eric.rasche.co.uk
#       for
#               Center for Phage Technology
#
# PODNAME: pause_histo.pl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Bio::DB::Sam;
use Bio::SeqIO;
use Data::Dumper;
use Statistics::Descriptive;

=head1 DESCRIPTION

PAUSE:histo is a tool to plot a complete historgram of number of aligned reads for every single base.

=head1 PARAMETERS

=over 4

=item C<bam_sense> bam file of reads mapping to sense strand

=item C<bam_antisense> bam file of reads mapping to antisense strand

=item C<genome> fasta sequence of the genome

=back

=cut

my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
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
		[
			'genome',
			'Genomic DNA as Fasta File',
			{ required => 1, validate => 'File/Input' }
		],
	],
	'outputs' => [
		[
			'pause_results',
			'Output TXT File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'pause_results.csv',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PAUSE_Histogram',
		'appname' => 'PAUSE:Histogram',
		'appdesc' =>
		  'plots a simple histogram of a bam mapping to a genome',
		'appvers' => '1.93',
	],
	'tests' => [
		#{
			#test_name    => "Default",
			#params => {
				#'sense_bam' => 't/Angus.sense.sam.bam.sorted.bam',
				#'antisense_bam' => 't/Angus.antisense.sam.bam.sorted.bam',
				#'genome' => 't/Angus.fa',
				#'pause_results_format' => 'YAML',
			#},
			#outputs      => {
				#'pause_results' => ['pause_results.yml','t/pause_results.yml' ],
			#}
		#},
	]

);

use CPT::Analysis::PAUSE;
my $pause = CPT::Analysis::PAUSE->new();

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
	name => 'pause_results',
	GGO => $ggo,
);

my $seqio_object = Bio::SeqIO->new( -file => $options->{genome}, -format => 'fasta');
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

	my %results;
	%results = $pause->histogram(
		data  => $sense_psam->coverage_density(),
		title => 'sense',
	);
	$output->varCRR(
		filename => sprintf( "%s.%s", $fasta_id, 'sense' ),
		extension => 'csv',
		data        => \%results,
		data_format => 'text/tabular',
		format_as   => 'CSV',
	);

	%results = $pause->histogram(
		data  => $antisense_psam->coverage_density(),
		title => 'antisense',
	);

	$output->varCRR(
		filename => sprintf( "%s.%s", $fasta_id, 'antisense' ),
		extension => 'csv',
		data        => \%results,
		data_format => 'text/tabular',
		format_as   => 'CSV',
	);
}
