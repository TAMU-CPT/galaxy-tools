#!/usr/bin/perl
#
#	   Code written by Eric Rasche
#		   mailto:rasche.eric@yandex.ru
#		   tel:   404.692.2048
#		   http://eric.rasche.co.uk
#	   for
#		   Center for Phage Technology
#

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my %non_hypothetical = ();

# The s@ implies an array of strings (i.e., the --file argument will apper multiple times)
my $options = $ggo->getOptions(
	'options' => [
		[
			'file|f',
			'Input file',
			{
				multiple => 1,
				required => 1,
				validate => 'File/Input',
			}
		],
		[
			'label',
			'Sheet name which corresponds to a given file',
			{ required => 1, multiple => 1, validate => 'String' }
		],
		[
			'genome|g', 'Input Genome
			File', { required => 1, validate => 'File/Input' }
		],
		[],
		[
			'dice_cutoff',
			'Dice',
			{
				required => 1,
				validate => 'Int',
				default  => 30,
			}
		],
		[
			'evalue_cutoff', 'Evalue',
			{ required => 1, validate => 'Int', default => -5 }
		],
		[
			'lookahead',
			'Lookahead for running ShineFind',
			{ required => 1, validate => 'Int', default => 30 }
		],
		[
			'disable_auto_hypothetical_novels|dahs',
"Disable the automatic labelling of hypothetical novels from Blast XML files based on the listing of Hits and Genes. If your FASTA header lines do not match your `/gene` tag in your `/products`, then you should set this flag"
		],
		[ 'full', 'Use Full output instead of Light output' ],
	],
	'outputs' => [
		[
			'annotation_table',
			'Annotation Table',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'AnnotationTable',
				data_format    => 'text/tabular',
				default_format => 'XLS'
			}
		],
	],
	'defaults' => [
		'appvers' => '1.94',
		'appid'   => 'AT_Table',
		'appname' => 'Annotation Table Generator',
		'appdesc' =>
'from a set of input BlastXML results, automatically filling in Hypothetical Novels',
	],
	'tests' => [
		{
			test_name => "Default",
			params => {
				'genome' => 'test-data/inputs/single.gbk',
				'file' => 'test-data/inputs/blast_results/blastp_bct.csv',
				'label' => 'BCT',
				'annotation_table_format' => 'YAML',
			},
			outputs => {
				'annotation_table' => ['AnnotationTable.yml', 'test-data/outputs/at_table2.deafult.yaml' ],
			}
		},
	],
);
run();

sub process_sheet_name {
	my ( $name, $i ) = @_;
	$name =~ s/[^A-Za-z0-9_]//g;

	# Just in case.
	if ( length($name) == 0 ) {
		$name = sprintf( 'sheet_%s', $i + 1 );
	}
	return $name;
}

sub run {
	print Dumper $options;
	my %results;

	my @files      = @{ $options->{file} };
	my @sheet_name = @{ $options->{label} };
	die 'must provide as many labels as files'
	  if scalar @files != scalar @sheet_name;

	# Loop across files
	for ( my $i = 0 ; $i < scalar(@files) ; $i++ ) {
		my $sheet_data_ref = process_file( $files[$i] );
		my $name = process_sheet_name( $sheet_name[$i], $i );
		$results{$name} = $sheet_data_ref;
	}
	$results{'Annotation Table'} = create_last_page();
	use CPT::OutputFiles;
	my $at_output = CPT::OutputFiles->new(
		name => 'annotation_table',
		GGO => $ggo,
	);
	$at_output->CRR(data => \%results);
}

sub process_file {
	my $response_ref = shift;
	print STDERR "Working on $response_ref\n";

	my @header_full = (
		"IDX",            #to sort with
		"Query ID",       #query id
		"Hit ID",         #NCBI link
		"NCBI Link",      # URL to ncbi   from ParentAccession
		"Description",    #Portal protein
		"Source",         #alldescs parsed
		"Hit Length",     #Lenght of hit
		"Bit Score",      # Bit score
		"evalue",
		"% dice",

		"query-start",
		"query-end",
		"subject-start",
		"subject-end",
		"identity",
		"positive",
		"gaps",

		"qseq,hseq"
	);

	my @data;
	my $orig_order = 1;

	# Fix 10e-50.0 to 10e-50
	if ( $options->{evalue_cutoff} =~ /\.0/ ) {
		$options->{evalue_cutoff} =~ s/\.0//g;
	}
	my $evalue_cutoff = "10e" . $options->{evalue_cutoff};

	open( my $csv_file, '<', $response_ref );
	while (<$csv_file>) {
		chomp $_;
		my (
			$qseqid,    $sseqid,  $pident, $length,
			$mismatch,  $gapopen, $qstart, $qend,
			$sstart,    $send,    $evalue, $bitscore,
			$sallseqid, $score,   $nident, $positive,
			$gaps,      $ppos,    $qframe, $sframe,
			$qseq,      $sseq,    $qlen,   $slen,
			$alldescs
		) = split( "\t", $_ );

		my $Identity = $pident * $qlen;
		my $dice = ( 2 * $Identity ) / ( $slen + $qlen );

		if (       $dice > $options->{dice_cutoff}
			&& $evalue <= $evalue_cutoff )
		{

			# Update non-hypo listing
			if ( !$options->{disable_auto_hypothetical_novels} ) {
				$non_hypothetical{$qseqid}++;
			}

			my @sseqid_arr = split( ';;', $sallseqid );
			my @gi_arr = map {
				if   ( $_ =~ /gi\|([0-9]+)\|/ ) { $1 }
				else                            { }
			} @sseqid_arr;

			my @unparsed_descs_arr = split( ';;', $alldescs );

			# transcriptional regulator [Paenibacillus alvei]
			my @descriptions = map {
				if   ( $_ =~ /^([^\[]*)/ ) { $1 }
				else                       { }
			} @unparsed_descs_arr;
			my @sources = map {
				if   ( $_ =~ /\[([^\]]*)\]/ ) { $1 }
				else                          { }
			} @unparsed_descs_arr;

#print STDERR $alldescs . "\t" . join("|",@descriptions) . "\t" . join("|",@sources) ."\n";

			my @template_row = (
				$qseqid,    # Query ID: Phagey200
				$sseqid
				,    # Hit ID: gi|491703197|ref|WP_005552565.1|
				'http://www.ncbi.nlm.nih.gov/protein/'
				,                         # NCBI Link
				'Portal Protein',         # description
				'Bacillus Phage Blah',    # Source
				$length,
				$bitscore,
				$evalue,
				$dice,
				$qstart,
				$qend,
				$sstart,
				$send,
				$Identity,
				$positive,
				$gaps,
				$qseq . $sseq,
			);
			for ( my $i = 0 ; $i < scalar(@descriptions) ; $i++ ) {
				$template_row[2] =
				  'http://www.ncbi.nlm.nih.gov/protein/'
				  . $gi_arr[$i];
				$template_row[3] = $descriptions[$i];
				$template_row[4] = $sources[$i];
				push( @data, [ $orig_order++, @template_row ] );
			}
		}
	}

	my %page_results = (
		'header' => \@header_full,
		'data'   => \@data,
	);
	return \%page_results;
}

sub create_last_page {
	my @data   = ();
	my @header = (
		'Gene',           'Strand',
		'Start Codon',    'Putative RBS',
		'Length',         'Putative Function',
		'Accession',      'Description',
		'Source',         'E-value',
		'Notes/comments', 'TMHMM',
		'LipoP',          'SignalP',
		'Protein Sequence'
	);

	# Use the naïve prediction algorithm that I use
	use CPT::Bio::RBS;
	my $rbs_predictor = CPT::Bio::RBS->new();
	$rbs_predictor->set_algorithm('naive');
	require Bio::SeqIO;
	require CPT::Bio;
	my $bio = CPT::Bio->new();

	my $seq = $bio->getSeqIO($options->{genome});
	my $seq_obj = $seq->next_seq;
	foreach my $feat_object ( $seq_obj->get_SeqFeatures ) {
		if ( $feat_object->primary_tag eq 'CDS' ) {
			my $name = $bio->_getIdentifier($feat_object);
			my $upstream = $bio->intelligent_get_seq($feat_object,
				upstream => 15,
				parent => $seq_obj,
			);
			my $upstream_for_analysis = lc(substr($upstream, 0, 10));
			my $upstream_missing = lc(substr($upstream, 10, 5));
			my @SDs = $rbs_predictor->predict($upstream_for_analysis);

			my @tmp;
			if ( $feat_object->strand == "1" ) {
				@tmp = (
					$name, '+',
					join("\n", map { $_->upstream . $upstream_missing } @SDs),
					join("\n", map { $_->rbs_seq } @SDs),
					abs(
						$feat_object->start -
						  $feat_object->end
					),
					(
						$non_hypothetical{$name}
						? ''
						: 'Hypothetical Novel'
					),    #'',# Assigned Putative Function
					      #'',# Assigned Acecssion
					      #'',# Assigned Desc
					      #'',# Assigned Evalue
					      #'',# Assigned Notes
					      #'',# Assigned TMHMM
				);
			}
			else {
				@tmp = (
					$name, '-',
					join("\n", map { $_->upstream . $upstream_missing } @SDs),
					join("\n", map { $_->rbs_seq } @SDs),
					abs(
						$feat_object->start -
						  $feat_object->end
					),
					(
						$non_hypothetical{$name}
						? ''
						: 'Hypothetical Novel'
					),
				);
			}
			push( @data, \@tmp );

			#push ($results{'Sheet1'}{'data'},\@tmp);
		}
	}

	my %results = (
		'header' => \@header,
		'data'   => \@data,
	);
	return \%results;
}

sub rc {
	my $val = shift;
	$val = reverse($val);
	$val =~ tr/actg/qzac/;
	$val =~ tr/qz/tg/;
	return $val;
}

=head1 NAME

Annotation Table Generator

=head1 DESCRIPTION

This script merges a number of CSV formatted blast hits, and produces an "annotation table" from that data. This annotation table consists of human-readable tabular blast results and the core annotation sheet which is used by students to record their putative gene assignments based on blast (and other) data.

=cut
