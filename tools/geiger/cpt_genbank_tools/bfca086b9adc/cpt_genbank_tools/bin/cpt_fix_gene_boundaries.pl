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
use Data::Dumper;

my $libCPT = CPT->new();

my $options = $libCPT->getOptions(
	'options' => [
		[
			'file|f',
			'Input file',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
	],
	'outputs' => [
		[
			'genbank',
			'GBK with Fixed Boundaries',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'boundaries',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank'
			}
		],
	],
	'defaults' => [
		'appid'   => 'FixGeneBoundaries',
		'appname' => 'Fix Gene Boundaries',
		'appvers' => '1.94',
		'appdesc' => 'to include RBS based on locus_tag',
	],
	'tests' => [
	],
);

my %args = ( 'file' => $options->file, );
use CPT::Bio;
my $bio = CPT::Bio->new();
my $seqobj = ${ $bio->requestCopy(%args) };

my %groups;
for my $feat ( $seqobj->get_SeqFeatures() ) {
	my $pt = $feat->primary_tag;
	if ( $pt eq 'RBS' ) {
		if ( $options->verbose ) {

 #printf("Looking at %s [%s]\n",$pt,$feat->get_tag_values('locus_tag'));
		}
		if ( $feat->has_tag('locus_tag') ) {
			my @locus_tag =
			  $feat->get_tag_values('locus_tag');
			my $preferred_locus_tag = $locus_tag[0];
			$groups{$preferred_locus_tag}{RBS} =
			  [ $feat->start, $feat->end ];
		}
		else {
			die 'Found an RBS without a locus tag!';
		}
	}
	elsif ( $pt eq 'gene' ) {
		if ( $options->verbose ) {

 #printf("Looking at %s [%s]\n",$pt,$feat->get_tag_values('locus_tag'));
		}
		if ( $feat->has_tag('locus_tag') ) {
			my @locus_tag =
			  $feat->get_tag_values('locus_tag');
			my $preferred_locus_tag = $locus_tag[0];
			$groups{$preferred_locus_tag}{gene} = $feat;
		}
		else {

			#print Dumper $feat;
			die 'Found an gene without a locus tag!';
		}
	}
}
foreach ( keys %groups ) {
	if ( defined $groups{$_}{RBS} ) {
		my ( $original_start, $original_end ) =
		  @{ $groups{$_}{RBS} };
		if ( defined $groups{$_}{gene} ) {
			if ( $options->{verbose} ) {
				printf STDERR "Updating start from %d to %d for %s\n", $original_start, $groups{$_}{gene}->start(),$_;
			}
			if ( $groups{$_}{gene}->strand() eq 1 ) {
				$groups{$_}{gene}->start($original_start);
			}
			else {
				$groups{$_}{gene}->end($original_end);
			}
		}
		else {
			warn "Could not find a gene feature associated with $_";
		}
	}
	else {
		warn "Issue with $_";
	}
}
use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
	name => 'genbank',
	libCPT => $libCPT,
);
$output->CRR(data => $seqobj);

=head1 DESCRIPTION

Some organisations expect gene boundaries to include RBSs and CDSs. In order to accomodate them, this tool modifies the boundaries of C<gene> features to include both relevant subfeatures. This is accomplished through use of the locus tag, so it MUST be identical between CDS, RBS, and gene.

=cut

