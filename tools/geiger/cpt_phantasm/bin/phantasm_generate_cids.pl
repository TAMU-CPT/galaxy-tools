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

# PODNAME: phantasm_generate_cids.pl
use CPT::GalaxyGetOpt;
use CPT::Util;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();
my $util = CPT::Util->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'file',
			'Input file',
			{
				required => 1,
				validate => 'File/Input'
			}
		],
	],
	'outputs' => [
		[
			'cassettes',
			'Identified gene cassettes',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'cassettes',
				data_format    => 'text/tabular',
				default_format => 'TSV_U',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PHAnTASM.CIDGeneration',
		'appname' => 'Cassette ID Generator',
		'appdesc' => 'calculate cassette ID string for genomes',
		'appvers' => '1.94',
	],
);

use File::ShareDir;
use File::Spec;
my $dir          = File::ShareDir::dist_dir('CPT-PHAnTASM');
my $cid_scheme   = File::Spec->catfile( $dir, 'clustering.yaml' );

my $seqin  = Bio::SeqIO->new( -file   => $options->{file} );
my %cassette_model = %{
	$util->JSONYAMLopts(
		'file'         => $cid_scheme,
	)
};

my %regex_containers;
my %custom_regexes;

foreach my $key ( keys %cassette_model) {
	foreach ( @{ $cassette_model{$key}{'members'} } ) {
		my $regex = qr/$_/i;
		$regex_containers{$regex} = $cassette_model{$key}{id};
	}
	foreach ( keys %{ $cassette_model{$key}{'custom'} } ) {
		$custom_regexes{$_} = {
			'title' => $cassette_model{$key}{id},
			'is'    => $cassette_model{$key}{'custom'}{$_}{'is'},
			'isnot' => $cassette_model{$key}{'custom'}{$_}{'isnot'},
		};
	}
}

use CPT::Bio;
my $bio = CPT::Bio->new();
my $seqio = Bio::SeqIO->new( -file => $options->{file}, -format => 'genbank' );

my $i = 0;

my @table_data;
while(my $seqobj = $seqio->next_seq()){
	my @cassette_model;

	foreach my $feat ( $seqobj->get_SeqFeatures ) {
		if ( $feat->primary_tag eq 'CDS' ) {
			my $feature_result = $bio->_getFeatureTag( $feat, 'product' ). $bio->_getFeatureTag( $feat, 'note' );
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
					my $ok_to_overwrite = 1;
					# Fixes a strange bug. If we have a
					# match, and then we have a subpart of
					# that match which hits a custom
					# element, we want to make sure that
					# it'll ONLY overwrite if we didn't
					# specifically exclude items like the
					# one we hit.
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
						$matched_result = $custom_regexes{$custom}{'id'};
					}
				}
			}
			if ($matched_result) {
				push(@cassette_model, [ $feat->strand == 1 ? '+': '-', $matched_result ]);
			}
		}
	}
	if(scalar(@cassette_model) > 4){
		my $raw_id = join('', map { join('',@{$_}) } @cassette_model);
		my ($orig, $mod) = collapse($raw_id);
		if(defined $mod){
			push(@table_data,[$seqobj->display_id(), $mod]);
		}
	}
}

my %data = (
	'Sheet1' => {
		header => ['#Genome', 'CID'],
		data => \@table_data,
	}
);

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
	name => 'cassettes',
	GGO => $ggo,
);
$output->CRR(data => \%data);



# Importing collapse.pl
sub collapse {
	my ($a) = @_;
	my @s = unpack("(A2)*", $a);
	if(scalar @s < 4){
		return;
	}
	my @operons = operons(@s);
	my @reduced;
	foreach(@operons){
		push(@reduced, algo(@$_));
	}
	my $mod = join('', map{join '', @$_} @reduced);
	my $orig = join('', map{join '', @$_} @operons);
	return ($orig, $mod);
}

sub operons {
	my @operons;
	my $cur_sign = '';
	my @c;
	foreach(@_){
		if(substr($_,0,1) eq $cur_sign){
			push(@c, $_);
		}else{
			if(scalar @c){
				my @copy;
				foreach(@c){push(@copy, "$_");}
				push(@operons,\@copy);
				@c = ();
			}
			$cur_sign = substr($_,0,1);
			push(@c, $_);
		}
	}
	my @copy;
	foreach(@c){push(@copy, "$_");}
	push(@operons,\@copy);
	return @operons;
}

sub algo {
	#my @d = window(2,2,@_);
	#my @e = rmrep(@d);
	my @e = rmrep(@_);

	return [@e];
}

sub rmrep {
	my @f;
	my $last = "";
	foreach(@_){
		if($_ ne $last){
			$last = $_;
			push(@f, $_);
		}
	}
	return @f;
}

=head2 window

Runs popcon for a window specified forwards/backwards.

=cut

sub window {
	my ($f,$r,@r) = @_;
	#print "window ($f,$r)> ". join("",@r) . "\n";
	my @e;
	my $called = 0;
	for($i=$r;$i<scalar @r - $f;$i++){
		my @z = @r[$i-$r .. $i+$f];
		push(@e, popcon(@z));
		$called = 1;
	}
	if(!$called){
		return @r;
	}
	return @e;
}

=head2 popcon

Incredibly simple algo which just chooses the most popular letter (named after
package popularity contest in linux/ubuntu)

=cut

sub popcon {
	my %h;
	#print 'popcon> ' .join('',@_), "\t";
	foreach(@_){
		$h{substr($_,1,1)}++;
	}
	#print Dumper \%h;
	my @top = sort { $h{$b} <=> $h{$a} } keys(%h);
	return substr($_[0],0,1) . $top[0];
}


=head1 NAME

PHAnTASM CID Generation Tool

=head1 DESCRIPTION

For an input B<multi-genome genbank> file, this tool will generate CIDs for every genome in the file.

What is a CID you ask? Why it's a genomic "cassette ID". We identify every cassette within the genome using a list of keywords that likely appear in the "note" or "product" tags, and generate a string which identifies this genome in a non-unique manner.

For a given genome with a replication, lysis, and morphogenesis cassette all on the same strand, you would receive a CID of "r+l+m+".

CIDS consist of a set of letters and +/- signs to denote strands. These strings can be used to compare the genomic layout of two phages.

=cut
