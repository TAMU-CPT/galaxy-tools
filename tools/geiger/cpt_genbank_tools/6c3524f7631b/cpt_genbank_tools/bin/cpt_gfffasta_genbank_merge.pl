#!/usr/bin/perl
#
#       Code written by Eric Rasche
#               mailto:rasche.eric@yandex.ru
#               tel:404.692.2048
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
			'file_fasta' => 'Fasta file that was used in the GFF producing tool',
			{
				required    => 1,
				validate    => 'File/Input',
				file_format => ['fasta'],
			}
		],
		[
			'file_gff' => 'GFF Output of tool',
			{
				required    => 1,
				validate    => 'File/Input',
			}
		],
		[
			'file_gbk' => 'Genbank file from which sequences were exported originally',
			{
				required    => 1,
				validate    => 'File/Input',
			}
		],
		['flatten_annotations', 'For hits which match to a feature, should their annotations be merged down into the parent feature\'s annotation? If this option is NOT set, then new features will be created for each of the hits' ],
	],
	'outputs' => [
		[
			'genbank',
			'Merged Genbank File',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'merged',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank'
			}
		],
	],
	'defaults' => [
		'appid'   => 'GFFGenbankMere',
		'appname' => 'Merge GFF output into Genbank',
		'appvers' => '1.94',
		'appdesc' =>
		  'given the gff output of a tool, the fasta input into that tool, and the original genbank file, this tool merges all three into a final genbank file'
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file_gff' => 'test-data/inputs/example.gff3',
				'file_fasta' => 'test-data/inputs/inputfa_for_gff3.fa',
				'file_gbk' => 'test-data/inputs/multi.gbk',
				'flatten_annotations' => '',
			},
			outputs => {
				'genbank' => ['merged.gbk', 'test-data/outputs/gff_fasta_gbk_merge.gbk'],
			}
		},
	],
);
# We use crc64 to match sequences.
use CPT::Util::CRC64;
my $crc = CPT::Util::CRC64->new();

# Load the fasta
use Bio::SeqIO;
open( my $fasta, '<', $options->{file_fasta} );
my $fasta_seqio = Bio::SeqIO->new( -fh => $fasta, -format => 'fasta' );
my ( $seq_id, $seq );
my %fasta_seq_map;
while ( my $seqobj = $fasta_seqio->next_seq() ) {
	$seq_id = $seqobj->display_id();
	$seq    = $seqobj->seq();

	my $hash = $crc->crc64($seq);
	# Can we assume unique sequences?
	$fasta_seq_map{$seq_id} = $hash;
}
close($fasta);

# Load/open GFF
use Bio::Tools::GFF;
open( my $fh, '<', $options->{file_gff} );
my $gffio = Bio::Tools::GFF->new( -fh => $fh, gff_version => 3);
my $feature;

# Associate the GFF features with fasta features and CRC hashes.
my %feature_map;
while ( $feature = $gffio->next_feature() ) {
	#$new_seq->add_SeqFeature($feature);
	my $hash = $fasta_seq_map{$feature->seq_id()};
	if(defined $hash){
		push(@{$feature_map{$hash}}, $feature);
	}else{
		warn "ERROR, could not find associated fasta for $feature";
	}
}
$gffio->close();

# Open Seqobj
use CPT::Bio;
my $bio = CPT::Bio->new();
use CPT::Bio::GFF_Parsing;
my $gff_conv = CPT::Bio::GFF_Parsing->new();
# Convert IPR121232 to Interpro:IPR123212
use CPT::Bio::Dbxref;
my $dbxref = CPT::Bio::Dbxref->new();

my $seqio_object = Bio::SeqIO->new(-file => $options->{file_gbk},-format=>'genbank');
my @updated_seqobjs;
while(my $seqobj = $seqio_object->next_seq()){
	foreach my $feat($seqobj->get_SeqFeatures()){
		my $seq = $bio->translate($bio->intelligent_get_seq($feat));
		# Remove stops
		$seq =~ s/[+*#]*//g;
		my $hash = $crc->crc64($seq);
		if(defined $feature_map{$hash}){
			# We have one or more features associated with this sequence
			my @matches = @{$feature_map{$hash}};
			my $strand = $feat->strand();
			my $start = $feat->start();
			my $end = $feat->end();
			foreach(@matches){
				if($options->{flatten_annotations}){
					# Copy all tags and values from sub-hits over.
					foreach my $tag($_->get_all_tags()){
						foreach my $value($_->get_tag_values($tag)){
							my $fixed = $gff_conv->fix_gff_tag($tag);
							if($fixed eq 'db_xref'){
								my @vals = $dbxref->get_prefix($value);
								if(scalar @vals == 1){
									$feat->add_tag_value($fixed, $vals[0] . ':' . $value);
								}elsif(scalar @vals > 1){
									warn "Found more than one possible match for $value: "
										. join(', ',@vals)
										. "\nIf any of those look 'more correct', please inform Eric\n";
									$feat->add_tag_value($fixed, $vals[0] . ':' . $value);
								}else{
									$feat->add_tag_value($fixed, $value);
								}
							}else{
								$feat->add_tag_value($fixed, $value);
							}
						}
					}
				}else{
					warn "Unimplemented"
				}
			}
		}
	}
	push(@updated_seqobjs, $seqobj);
}

# Loop over features

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
	name => 'genbank',
	libCPT => $libCPT,
);
$output->CRR(data => \@updated_seqobjs);

=head1 DESCRIPTION


=cut
