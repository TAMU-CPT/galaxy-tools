#!/usr/bin/perl
#
#	   Code written by Eric Rasche
#		  mailto:rasche.eric@yandex.ru
#		  tel:   404.692.2048
#		  http://eric.rasche.co.uk
#	   for
#		  Center for Phage Technology
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
			'Input File',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
		[
			'clustered_by_tag',
'If items in the genome are already clustered with a specific tag'
			  . '(i.e.; gene, CDS, RBS all share a locus_tag), then we can use this information to help renumber',
			{ validate => 'String' }
		],
		[
			'tag_to_update', 'Which tag is used to store gene
			numbering',
			{
				validate => 'String',
				default  => 'locus_tag'
			}
		],
		[
			'string_prefix', 'A string to use as
		a prefix for the numbering. Will be used as XXXXXNNN' . 'where
		XXXXX is the string. Using "display_id" has special meaning: it
		will use the genome\'s' . ' name/accession number.',
			{
				validate => 'String',
				default  => 'display_id'
			}
		],
		[
			'leading_zeros',
'Number of leading zeros. The number will be padded to this length',
			{ validate => 'Int', default => 3 }
		],
	],
	'outputs' => [
		[
			'genbank',
			'Renumbered Genbank File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'renumbered',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank',
			}
		],
		[
			'change_table',
			'Table of gene name changes',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'renumbered',
				data_format    => 'text/tabular',
				default_format => 'CSV',
			}
		],
	],
	'defaults' => [
		'appid'   => 'RelabelGBKTags',
		'appname' => 'Relabel Genbank Tags',
		'appvers' => '1.94',
		'appdesc' => 're-labels genbank tags according to rules',
	],
	'tests' => [
		{
			test_name => "Default",
			params => {
				'file' => 'test-data/inputs/single.gbk',
				'string_prefix' => 'gp',
			},
			outputs => {
				'genbank' => ['renumbered.gbk', 'test-data/outputs/renumbered.gbk'],
				'change_table' => ['renumbered.Sheet1.csv', 'test-data/outputs/renumbered.csv'],
			}
		},
		{
			test_name => "cluster by tag 'locus_tag'",
			params => {
				'file' => 'test-data/inputs/single.gbk',
				'string_prefix' => 'gp',
				'clustered_by_tag' => 'locus_tag',
			},
			outputs => {
				'genbank' => ['renumbered.gbk', 'test-data/outputs/renumbered.cluster.gbk'],
				'change_table' => ['renumbered.Sheet1.csv', 'test-data/outputs/renumbered.cluster.csv'],
			}
		},
	],
);

my @clusters;

#sub which_tag {
#	my %tag_usage_map;
#	my %reverse_tag_map;
#	# (
#	#   gene3 => {
#	#	   tags => {
#	#		   locus_tag = "gene009",
#	#		   Note = "blah lahblbabklh",
#	#	   },
#	#	   reverse_tags {
#	foreach my $gene(@genes){
#		# add gene to tag usage map (start_stop_strand)
#		# calculate avg distance for each tag
#	}
#}
#
#
#use List::Util qw/min/;
#sub LevenshteinDistance{
#	my ($s,$len_s,$t,$len_t) = @_;
#	if(length($len_s) == 0){ return $len_t; }
#	if(length($len_t) == 0){ return $len_s; }
#
#	my $cost;
#	if(substr($s,$len_s-1,1) eq substr($t,$len_t-1,1)){
#		$cost = 0;
#	}else{
#		$cost = 1;
#	}
#
#	my $a = LevenshteinDistance($s,$len_s-1,$t,$len_t) + 1;
#	my $b = LevenshteinDistance($s,$len_s, $t,$len_t -1 ) + 1;
#	my $c = LevenshteinDistance($s,$len_s - 1,$t,$len_t -1) + $cost;
#	return min($a,$b,$c);
#}

use CPT::Bio;
my $bio = CPT::Bio->new();
my $seqobj = ${ $bio->requestCopy( file => $options->{file} ) };

# Global, yuck!
my $idx = 0;

my $string_prefix;
if ( $options->{string_prefix} && $options->{string_prefix} ne 'display_id' ) {
	$string_prefix = $options->{string_prefix};
}
else {
	$string_prefix = $seqobj->display_id();
}
my $format_string = $string_prefix . '%0' . $options->{leading_zeros} . 'd';

# If we're using clustering by tag (i.e., the user has specified a tag)
#if($options->{clustered_by_tag}){
#my %feat_clust;
#foreach my $feat ($seqobj->get_SeqFeatures){
#if($feat->has_tag($options->{clustered_by_tag})){
#my @val = $feat->get_tag_values($options->{clustered_by_tag});
#push(@{$feat_clust{$val[0]}}, $feat);
#}
#}
#foreach(keys %feat_clust){
#my @feats = @{$_};
#foreach my $feat(@feats){
#update_feature($feat);
#}
#}
#}else{
## Otherwise we just go through sequentially, relying on BioPerl for correct data. :/
#foreach my $feat ( $seqobj->get_SeqFeatures){
#update_feature($feat);
#}
#}

my %feature_cluster;

sub attempt_cluster {
	my ($feat) = @_;

# If this feature is already known to the cluster by its tag...
#printf "Attempting to cluster [%s %s %s %s]\n" , $feat->display_name(), $feat->start(), $feat->end(), $feat->strand();
	if ( $feat->has_tag( $options->{clustered_by_tag} ) ) {

		#printf "\thas tag\n";
		my @val = $feat->get_tag_values( $options->{clustered_by_tag} );
		if ( $feature_cluster{ $val[0] } ) {

	      # We have a hit, should probably just cluster it(?)
	      # Better yet would be to double check by looking for overlap > 80%
			push( @{ $feature_cluster{ $val[0] } }, $feat );
			return;
		}
		else {
			$feature_cluster{ $val[0] } = [$feat];
		}
	}

	# Loop through existing clusters
	foreach my $gene_arrays ( keys(%feature_cluster) ) {
		my @gene_set = @{ $feature_cluster{$gene_arrays} };
		foreach my $clustered_feat (@gene_set) {
			if (       $feat->start() == $clustered_feat->start()
				|| $feat->end() == $clustered_feat->end() )
			{
				push( @gene_set, $feat );
				$feature_cluster{$gene_arrays} = \@gene_set;
				return;
			}
		}
	}

	# Otherwise, no existing tags, nor matching start/ends with anything
	# already in there
	push( @{ $feature_cluster{ sprintf( 'cluster_%s', rand() ) } }, $feat );

	# We don't really care how we identify it, as we're going to do the
	# renaming by start/end anyway
}

# Here we update the keys to be the lowest value in the array (i.e., the left
# most gene start/end). Seems like a reasonable way to organise that stuff.
sub add_start_info {
	my $lidx = 0;
	foreach my $key ( keys(%feature_cluster) ) {

		# Set to zero, will reset immediately
		my $left_most = 0;

		# Looping across features in this group
		foreach ( @{ $feature_cluster{$key} } ) {
			if ( $left_most == 0 ) {
				$left_most = $_->start;
			}

			# Update if it's less than
			if ( $_->start < $left_most ) {
				$left_most = $_->start;
			}
			if ( $_->end < $left_most ) {
				$left_most = $_->end;
			}
		}
		if ( $left_most != 0 ) {

#print "[$left_most] Adjusting name for $key to " .  sprintf('__%s__' , $lidx) . "\n";
			$feature_cluster{ sprintf( '__%s__', $left_most ) } =
			  $feature_cluster{$key};
			$lidx++;
		}
		else {
			warn 'THIS IS BAD';
		}
	}
}

# The actual clustering
foreach my $feat ( $seqobj->get_SeqFeatures ) {
    if($feat->primary_tag eq 'CDS' || $feat->primary_tag eq 'RBS' || $feat->primary_tag eq 'gene'){
	attempt_cluster($feat);

}}
add_start_info();

my @good_keyset = map {
	if   ( $_ =~ /__([0-9]+)__/ ) { $1; }
	else                          { }
} keys(%feature_cluster);
foreach my $cluster ( sort { $a <=> $b } @good_keyset ) {

     #print "Cluster $cluster => " . $feature_cluster{"__${cluster}__"} . " \n";
	renumber_cluster( @{ $feature_cluster{"__${cluster}__"} } );
}

sub renumber_cluster {
	my (@feats) = @_;
	foreach (@feats) {
		update_feature($_);
	}
	$idx++;
}

my @data;
sub update_feature {
	my $feat = shift;

	my @previous_values = (sprintf('%s %s %s', $feat->start(), $feat->end(), $feat->strand()));
	if ( $feat->has_tag( $options->{tag_to_update} ) ) {
		@previous_values = $feat->get_tag_values($options->{tag_to_update});
		$feat->remove_tag( $options->{tag_to_update} );
	}

	push(@data,[$previous_values[0],sprintf($format_string, $idx)]);
	$feat->add_tag_value( $options->{tag_to_update},
		sprintf( $format_string, $idx ) );
}

my %fd = (
	'Sheet1' => {
		header => ['Original Name','New Name'],
		data => \@data,
	}
);


use CPT::OutputFiles;
# Sequence
my $crr_output = CPT::OutputFiles->new(
	name => 'genbank',
	libCPT => $libCPT,
);
$crr_output->CRR(data => $seqobj);
# Delta Table
$crr_output = CPT::OutputFiles->new(
	name => 'change_table',
	libCPT => $libCPT,
);
$crr_output->CRR(data => \%fd);

=head1 DESCRIPTION

Tool attempts to modify GenBank files to relabel all tags in the genome according to user-specified format. It is pending a re-write which will fix several important bugs

=cut
