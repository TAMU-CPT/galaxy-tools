#!/usr/bin/perl
#
#       Code written by Eric Rasche
#               mailto:rasche.eric@yandex.ru
#               tel:   404.692.2048
#               http://eric.rasche.co.uk
#       for
#               Center for Phage Technology
#

use strict;
use warnings;

use CPT;
use CPT::Plot::ArtemisColours;
use Data::Dumper;

my $libCPT = CPT->new();
my %modes  = (
	'remove'  => 'Remove all other color tags before inserting new one',
	'append'  => 'Append new color tags, leaving the old in place',
	'prepend' => 'Prepend new color tags, leaving the old in place',
	'ignore'  => 'Do nothing, use the old color tag',
);

my $options = $libCPT->getOptions(
	'options' => [
		[
			'file',
			'Input file',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
		[],
		['Behavior'],
		[
			'mode',
			=> "How should we treat the addition of new color tags to genes with existing color tags?",
			{
				required => 1,
				validate => 'Option',
				options  => \%modes
			}

		],
	],
	'outputs' => [
		[
			'genbank',
			'Recolored Genbank File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'colored',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank',
			}
		],
	],
	'defaults' => [
		'mode'    => 'remove',
		'remove'  => 1,
		'appid'   => 'Color',
		'appname' => 'Color',
		'appdesc' =>
'Add a `/color` tag to GenBank files based on text in the `/product` field.',
		'appvers' => '1.94',
	],
	'tests' => [
	]
);

use CPT::Bio;
my $bio = CPT::Bio->new();
my $seqin  = $bio->getSeqIO($options->{file});
my $seqout = Bio::SeqIO->new( -format => 'genbank' );

use File::ShareDir;
use File::Spec;
my $data_dir = File::ShareDir::dist_dir('CPT-CLI');
my $color_model_location = File::Spec->catfile(
	$data_dir,
	'colouring.yaml'
);
# Check for their color_model, otherwise use a builtin(?)
my %color_model;
%color_model = %{
	$libCPT->JSONYAMLopts(
		'file'         => $color_model_location,
		'ext_override' => 'yaml'
	)
  };

my %regex_containers;
my %custom_regexes;

my $colour_parser = CPT::Plot::ArtemisColours->new( format => 'artemis' );

foreach my $key ( keys %color_model ) {
	foreach ( @{ $color_model{$key}{'members'} } ) {
		my $regex = qr/$_/i;
		$regex_containers{$regex} =
		  $colour_parser->getColour( $color_model{$key}{'color'} );
	}
	foreach ( keys %{ $color_model{$key}{'custom'} } ) {
		$custom_regexes{$_} = {
			'title' => $color_model{$key}{'title'},
			'color' => $colour_parser->getColour(
				$color_model{$key}{'color'}
			),
			'is'    => $color_model{$key}{'custom'}{$_}{'is'},
			'isnot' => $color_model{$key}{'custom'}{$_}{'isnot'},
		};
	}
}

# Only process the first sequence object. TODO
my $seqobj = ${ $bio->requestCopy( file => $options->{file} ) };
foreach my $feat ( $seqobj->get_SeqFeatures ) {
	if ( $feat->primary_tag eq 'CDS' ) {
		my $feature_result = $bio->_getFeatureTag( $feat, 'product' )
		. $bio->_getFeatureTag( $feat, 'note' );
		my $matched_result = '';

		foreach my $regex ( keys %regex_containers ) {
			if ( $feature_result =~ $regex ) {
				$matched_result = $regex_containers{$regex};
			}
		}
		foreach my $custom ( keys %custom_regexes ) {

			# Do we possibly want to have another look at this one?
			my $care = 0;
			foreach ( @{ $custom_regexes{$custom}{'is'} } ) {
				if ( $feature_result =~ /\b$_\b/i ) {
					$care = 1;
				}
			}
			foreach ( @{ $custom_regexes{$custom}{'isnot'} } ) {
				if ( $feature_result =~ /\b$_\b/i ) {
					$care = 0;
				}
			}
			if ($care) {
				my $is_okay         = 1;
				my $ok_to_overwrite = 1
				  ; # Fixes a strange bug. If we have a match, and then we have a subpart of that match which hits a custom element, we want to make sure that it'll ONLY overwrite if we didn't specifically exclude items like the one we hit.
				foreach ( @{ $custom_regexes{$custom}{'is'} } )
				{
					if ( $feature_result !~ /\b$_\b/i ) {
						$is_okay = 0;
					}
				}
				foreach (
					@{ $custom_regexes{$custom}{'isnot'} } )
				{
					if ( $feature_result =~ /\b$_\b/i ) {
						$ok_to_overwrite = 0;
						$is_okay         = 0;
					}
				}
				if ( $is_okay && $ok_to_overwrite ) {
					$matched_result =
					  $custom_regexes{$custom}{'color'};
				}
			}
		}
		my $matched_colour = $matched_result;
		if ($matched_colour) {

			#print "Matched $matched_colour\n";
			my @colors;
			if ( $options->{'mode'} eq 'remove' ) {

			    # Don't need to muck about, just set it as our list.
				push( @colors, $matched_colour );
			}
			elsif (    $options->{'mode'} eq 'append'
				|| $options->{'mode'} eq 'prepend' )
			{
				if ( $feat->has_tag('color') ) {
					@colors =
					  $feat->get_tag_values('color');
					if ( $options->{'mode'} eq 'append' ) {
						push( @colors,
							$matched_colour );
					}
					else {
						unshift( @colors,
							$matched_colour );
					}
				}
				else {
					@colors = ($matched_colour);
				}
			}
			elsif ( $options->{mode} eq 'ignore' ) {
			}
			if ( $feat->has_tag('color') ) {
				$feat->remove_tag('color');
			}
			foreach (@colors) {
				$feat->add_tag_value( 'color', $_ );
			}
		}
	}
}

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
        name => 'genbank',
        libCPT => $libCPT,
);
$output->CRR(data => $seqobj);

=head1 NAME

Genbank Coloring Utility

=head1 DESCRIPTION

In order to easily visualize clusters of genes in phages with related functions, the CPT has found it of use to be able to automatically apply colors to genes which can then be visualized in Artemis and the CPT's Genome Mapper. These coloring rules are based off of analysis of common tags for phage genomes.

=over 4



=item Biosynthesis

    Contains text: "Deoxynucleoside", "Deoxyribonucleotidase", "Deoxyuridine", "Ribonucleoside-diphoshate reductase", "Serine kinase", "threonine kinase", "cytidine deaminase", "dUTPase", "dUTPase", "deoxynucleotide", "dihydrofolate reductase", "glutaredoxin", "guanylate kinase", "reductase", "ribonucleotide reductase", "thioredoxin", "thymidylate"

Color: L<rgb(0,0,255)|https://duckduckgo.com/?q=rgb+0+0+255>

=item Defense

    Contains text: "rII", "rIIA", "rIIB", "rex", "rexA", "rexB", "ocr", "dar", "darA", "darB"

Color: L<rgb(200,255,200)|https://duckduckgo.com/?q=rgb+200+255+200>

=item DNA Packaging

    Contains text: "terminase"

Color: L<rgb(0,255,255)|https://duckduckgo.com/?q=rgb+0+255+255>

=item DNA Replication/Recombination

    Contains text: "Clamp", "DNA binding protein", "DNA end Protector", "DNA ligase", "DexA", "DnaA", "DnaB", "DnaQ", "Helicase", "RNA ligase", "RNaseH", "RecA", "RecF", "Recombination", "RuvC", "UvsW", "UvsY", "helicase", "holliday junction", "phosphoesterase", "primase", "recombinase", "recombination", "repair", "single strand annealing", "topoisomerase", "whisker", "sliding", "methylase", "methyltransferase", "mom", "glucosyl\s*transferase", "glycosyl\s*transferase", "integrase"

    Special rules:
        - "polymerase" but not "rna polymerase" or "polymerase sigma factor"
        - "nuclease" but not "HNH" or "homing endonuclease"

Color: L<rgb(255,255,0)|https://duckduckgo.com/?q=rgb+255+255+0>

=item HNH/Homing/GIY-YIG

    Contains text: "HNH", "homing endonuclease", "GIY-YIG"

Color: L<rgb(200,150,100)|https://duckduckgo.com/?q=rgb+200+150+100>

=item Lysis

    Contains text: "antiholin", "holin", "endolysin", "spanin", "peptidoglycan", "amidase", "transglycosylase", "carboxypeptidase"

    Special rules:
        - "lysozyme" but not "lysozyme baseplate" or "tail lysozyme"

Color: L<rgb(255,0,255)|https://duckduckgo.com/?q=rgb+255+0+255>

=item Morphogenesis

    Contains text: "tail\s*spike", "fiber", "neck", "sheath", "tube", "pectin", "prohead", "scaffold", "capsid", "head", "head-to-tail joining", "pre-neck", "Tape", "tailspike", "structural", "morphogenesis", "assembly", "chaperone", "joining", "decoration", "protease", "frameshift", "portal"

    Special rules:
        - "tail" but not "tail lysozyme"
        - "baseplate" but not "lysozyme baseplate"

Color: L<rgb(135,206,250)|https://duckduckgo.com/?q=rgb+135+206+250>

=item Novel

    Contains text: "Novel"

Color: L<rgb(170,170,170)|https://duckduckgo.com/?q=rgb+170+170+170>

=item Regulation

    Contains text: "Translational repressor", "RegA", "RegB", "regulatory", "regulator", "transcriptional repressor", "anti-repressor", "rna", "rna polymerase", "Sigma Factor"

Color: L<rgb(255,165,0)|https://duckduckgo.com/?q=rgb+255+165+0>

=item terminator

    Contains text: "terminator"

Color: L<rgb(0,255,0)|https://duckduckgo.com/?q=rgb+0+255+0>

=item tRNAs

    Contains text: "tRNA"

Color: L<rgb(113,188,120)|https://duckduckgo.com/?q=rgb+113+188+120>

=back

=cut
