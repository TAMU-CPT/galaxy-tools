#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo = CPT::GalaxyGetOpt->new();

my $options = $ggo->getOptions(
	'options' => [
		[
			'gff3',
			'Input GFF3 file',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['gff3'],
			},
		],
		[
			'operation',
			'Operation to do on GFF3 file',
			{
				required => 1,
				validate => 'Option',
				options => {
					'flip' => 'Flip strand of GFF features',
					'invert' => 'Change coding regions to non-coding regions, and vice versa',
				},
			}
		],
		[
			'length',
			'Genome length',
			{
				validate => 'Int',
			}
		],
	],
	'outputs' => [
		[
			'gff3',
			'Modified GFF3 File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'modified.gff3',
				data_format    => 'genomic/interval',
				default_format => 'GFF3'
			}
		],
	],
	'defaults' => [
		'appid'   => 'edu.tamu.cpt.gff.GffEditor',
		'appname' => 'GFF Editor',
		'appvers' => '1.94',
		'appdesc' => 'modify GFF3 files',
	],
	'tests' => [
	],
);
use Bio::Tools::GFF;

use CPT::Bio;
my $bio = CPT::Bio->new();

my $gff = new Bio::Tools::GFF( -file => $options->{gff3}, -gff_version => 3 );
my $gffio = Bio::Tools::GFF->new(-gff_version => 3);
my @features;

if($options->{operation} eq 'flip'){
	while ( my $feature = $gff->next_feature ) {
		if($feature->strand > 0){
			$feature->strand(-1)
		}
		else{
			$feature->strand(1)
		}
		push(@features, $feature);
	}
}
elsif($options->{operation} eq 'invert'){
	my %coverage;
	foreach(1..$options->{length}){
		$coverage{$_} = 0;
	}

	while ( my $feature = $gff->next_feature ) {
		my $loc = $feature->location;
		if(ref($loc) eq 'Bio::Location::Simple'){
			foreach(my $i = $feature->start; $i < $feature->end; $i++){
				$coverage{$i} = 1;
			}
		}else{
			for my $location ( $loc->sub_Location ) {
				foreach(my $i = $feature->start; $i < $feature->end; $i++){
					$coverage{$i} = 1;
				}
			}
		}
	}

	my $last_state = $coverage{1};
	my $string='';
	for(my $i = 1; $i<$options->{length};$i++){
		$string .= $coverage{$i};
	}

	while($string =~ /([0]+)/g){
		my $start = $-[0];
		my $end = $+[0];
		my $f = Bio::SeqFeature::Generic->new(
			-start => $start,
			-end => $end,
			-strand => 1,
			-primary_tag => 'CDS',
		);
		push(@features, $f);
	}
}

my $output = "##gff-version 3\n";
# Print features
foreach my $feature(@features){
	$output .= $feature->gff_string($gffio) . "\n";
}

use CPT::OutputFiles;
my $crr_output = CPT::OutputFiles->new(
	name => 'gff3',
	GGO => $ggo,
);
$crr_output->CRR(data => $output);

=head1 DESCRIPTION

Given a GFF3 set of annotations and a fasta file, this tool will create a single genbank file out of them

=cut
