#!/usr/bin/perl
#
#       Code written by Eric Rasche
#               mailto:rasche.eric@yandex.ru
#               tel:   404.692.2048
#               http://eric.rasche.co.uk
#       for
#               Center for Phage Technology
#

# PODNAME: cpt_genome_mapper.pl

use strict;
use warnings;

use CPT::Util;
use CPT::GalaxyGetOpt;
use CPT::Plot::Base;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();
my $cpt_util = CPT::Util->new();

my %justification = (
	justify =>
'Justify the strands of the genome, just like with text in a word processor',
	leftalign =>
'Left align the strands of the genome, the right edges will be uneven',
);
my %view = (
	"alt_every"  => 'Alternate each gene up/down, like DNA Master',
	"alt_random" => 'Alternate genes up/down randomly',
	"alt_none"   => 'No alternation. Genes are plotted on the same line',
	"alt_artemis" =>
	  'Alternate genes up/down by reading frame, like Artemis',
);

my $options = $ggo->getOptions(
	'options' => [
		[
			'file',
			'Input file',
			{
				required => 1,
				validate => 'File/Input',
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
		[],
		[
			'double_line_for_overlap',
'Use a double line where a section of the genome is plotted twice',
		],
		[
			'justification' =>
'How should each strand of the genome be treated when plotting',
			{
				required => 1,
				validate => 'Option',
				options  => \%justification,
				default  => 'justify'
			}
		],
		[
			'separate_strands',
'+/- strand should be plotted on different tracks (artemis style)',
		],
		[
			'view',
'This option controls the arrangement of genes in a track',
			{
				required => 1,
				validate => 'Option',
				options  => \%view,
				default  => 'alt_artemis'
			}
		],

		[
			'x_offset',
			'The border on the X-axes',
			{
				required => 1,
				validate => 'Int',
				min      => 0,
				default  => 30,
			}
		],
		[
			'y_offset',
			'The border on the Y-axes',
			{
				required => 1,
				validate => 'Int',
				min      => 0,
				default  => 100,
			}
		],
		[
			'opacity',
			'Opacity of the plotted genes',
			{
				required => 1,
				validate => 'Float',
				min      => 0,
				max      => 1,
				default  => 0.7,
			}
		],
		[],
		[
			'rows',
			'Number of rows the genome should be split into',
			{
				required => 1,
				validate => 'Int',
				min      => 1,
				default  => 3,
			}
		],
		[
			'split_factor',
'Factor by which each line overcompensates. If you have a short piece of the genome at the end of your last row, bump this up until each line of the genome is approximately equivalent.',
			{
				required => 1,
				validate => 'Float',
				min      => 1,
				max      => 1.7,
				default  => 1.02,
			}
		],
		[
			'width_mode',
			'How to size the plot horizontally.',
			{
				required => 1,
				validate => 'Option',
				options => {'dynamic' => 'Dynamically resize based on a zoom factor', 'static' => 'Static map width' },
				default  => 'dynamic',
			}
		],
		[
			'width_value',
			'Adjusts width of the plot using behaviour dependent on width_mode.',
			{
				required => 1,
				validate => 'Int',
				min      => 1,
				default  => 80,
			}
		],
		[
			'inter_line_separation|ils',
			'Separation between each consecutive plotted row',
			{
				required => 1,
				validate => 'Int',
				min      => 1,
				default  => 200,

			}
		],
		[],
		['Labelling'],
		[
			'label',
			'Enable Labelling',
			{
				validate => 'Flag',
				default  => 0,
			}
		],
		[
			'label_position',
			'Label position',
			{
				validate => 'Option',
				options  => {
					'above' =>
					  'Plot labels above genomic features',
					'on' =>
'Plot labels on top of genomic features'
				},
				default => 'above',
			}
		],
		[
			'label_shrink_mode',
			'How should the label be treated if it is too large?',
			{
				validate => 'Option',
				options  => {
					'shrink' =>
'Reduce font size to make it generally proportional to genomic feature size',
					'cutoff' =>
'Anything labels which would not fit are not printed'
				},
				default => 'cutoff',
			}
		],
		[
			'label_callouts',
			'Use callouts if labelling above'
		],
		[
			'label_from',
			'Where should label text be obtained from?',
			{
				validate => 'Option',
				options  => {
					'numeric' =>
'Basic numeric labels, which increase sequentially for selected features (must specify label_numeric_features)',
					'custom' =>
'Use a custom query (must specify label_query)',

				},
				default => 'numeric',
			}
		],
		[
			'label_text_source',
'Which tag (e.g., note, locus_tag, product) should be used as label text?',
			{
				validate => 'String',
				default  => 'locus_tag',
			}
		],
		[
			'label_numeric_features',
'If you selected to label with numbers, specify this option as many times as needed to select the features you wish to number',
			{
				validate => 'Genomic/Tag',
				multiple => 1,
			}
		],
		[
			'label_query',
'If you selected to label from a custom query, specify that here. An example query would be: !contains:"Hypothetical" key:"CDS,tRNA" tag:"locus_tag" which would search all CDSs and tRNAs which have a locus_tag and do not contain the text "Hypothetical"',
			{ validate => 'String', }
		],
	],
	'outputs' => [
		[
			'genome_map',
			'SVG Genome Map',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'map',
				data_format    => 'image/svg',
				default_format => 'SVG',
			}
		],
	],
	'defaults' => [
		'appid'   => 'Genome_Map',
		'appname' => 'Genome Mapper',
		'appdesc' =>
		  'Maps genomes into SVG files, similar to DNA Master',
	]
);

if(defined $options->{label_query}){
	$options->{label_query} =~ s/__dq__/"/g;
}
use Bio::SeqIO;
my $seqio = Bio::SeqIO->new(-file => $options->{file}, -format => 'Genbank');
my $seq = $seqio->next_seq();
my @features = $seq->get_SeqFeatures();
my $genome_length = $seq->length();

use File::ShareDir;
use File::Spec;
my $dir          = File::ShareDir::dist_dir('CPT-Annotation-Mapper');
#my $dir = '../data/';
my $color_scheme = File::Spec->catfile( $dir, 'cpt-color-scheme.json' );
my %color_scheme = %{ $cpt_util->JSONYAMLopts( 'file' => $color_scheme ) };

my %feat_group_ref = ();
my $maxRowLength;
my %rowdata     = ();
my $svg_control = CPT::Plot::Base->new(
	'double_line_for_overlap' => $options->{'double_line_for_overlap'},
	'justified'               => $options->{'justification'},
	'separate_strands'        => $options->{'separate_strands'},
	'view'                    => $options->{'view'},
	'x_offset'                => $options->{'x_offset'},
	'y_offset'                => $options->{'y_offset'},
	'opacity'                 => $options->{'opacity'},
	'rows'                    => $options->{'rows'},
	'width_mode' => $options->{'width_mode'},
	'width_value' => $options->{'width_value'},
	'ils' => ( $options->{'separate_strands'} ? 1 : 0.5 ) *
	  $options->{'inter_line_separation'},
	'query'        => $options->{'label_with_tags_custom'},
	'split_factor' => $options->{'split_factor'},

	'label'                  => $options->{'label'},
	'label_pos'              => $options->{'label_position'},
	'label_shrink_mode'      => $options->{'label_shrink_mode'},
	'label_callouts'         => $options->{'label_callouts'},
	'label_from'             => $options->{'label_from'},
	'label_text_source'      => $options->{'label_text_source'},
	'label_numeric_features' => $options->{'label_numeric_features'},
	'label_query'            => $options->{'label_query'},

	'color_scheme' => \%color_scheme,
	'rows'         => $options->{rows},

	'features'     => \@features,
	'genome_length'=> $genome_length,
	#'file'         => $options->{file},
	'_ft_count'              => 0,   #Used for alternating the items up/down
	'_internal_maxrowlength' => 0,
);
$svg_control->init();
$svg_control->partitionLines();
$svg_control->createSVG();

use CPT::OutputFiles;
my $map = CPT::OutputFiles->new(
	name => 'genome_map',
	GGO => $ggo,
);
$map->CRR(data => $svg_control->getSVG());
