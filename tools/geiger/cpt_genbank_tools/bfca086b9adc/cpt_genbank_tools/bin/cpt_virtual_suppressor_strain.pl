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
use Storable qw/dclone/;

my $libCPT  = CPT->new();
my $options = $libCPT->getOptions(
	'options' => [
		[ 'file', 'Input file',
			{
				validate => 'File/Input',
				required => 1,
			}
		],
		[ 'suppress', 'Suppress this stop codon. All features with this codon will be extended to the next available stop codon.',
			{
				validate => 'Option',
				required => 1,
				multiple => 1,
				options => {
					'amber' => 'TAG (Amber)',
					'ochre' => 'TAA (Ochre)',
					'opal'  => 'TGA (Opal/Umber)',
				}
			}
		],
		[ 'replacement', 'Replace a suppressed stop codon with a user specified amino acid. Defaults to an N (because of bioperl restrictions), but this should be replaced with an 1 letter amino acid. This option is ONLY used if you also specify to translate the sequence below',
			{
				validate => 'String',
				required => 1,
				multiple => 1,
				default => ['N'],
			}
		],
		[ 'translate', 'Translate to amino acids', { validate => 'Flag' }],
	],
	'outputs' => [
		[
			'results',
			'Genbank File with suppressed stops',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'suppressed',
				data_format    => 'genomic/annotated',
				default_format => 'Genbank',
			}
		],
		[
			'modified_features',
			'Fasta File with modified features',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'modified',
				data_format    => 'genomic/raw',
				default_format => 'Fasta',
			}
		],
	],
	'defaults' => [
		'appid'   => 'VirtSupStr',
		'appname' => 'Virtual Suppressor Strain',
		'appdesc' => 'generates a genome with specified stop codons suppressed',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default GBK",
			params => {
				'file' => 'test-data/inputs/multi.gbk',
				'suppress' => 'amber',
				'translate' => '',
			},
			outputs => {
				'results' => ['suppressed.gbk', 'test-data/outputs/virt_suppressor_strain.gbk'],
				'modified_features' => ['modified.fa', 'test-data/outputs/virt_suppressor_strain.fa'],
			},
		},
	],
);



# Iterate over lists to create map
use List::MoreUtils qw/each_array/;
my $it = each_array(@{$options->{suppress}}, @{$options->{replacement}});
my %tmp;
while(my ($sup, $rep) = $it->()){
	$tmp{$sup} = $rep;
}

# Make sure they haven't specified ALL the stops
my @stops = @{$options->{suppress}};
my @uniq_stops = keys(%tmp);
if(scalar(@uniq_stops) > 2){
	die 'Cannot specify all stops. Features will have nowhere to stop!';
}

my %stop_ignore_map = (
	'TAG' => defined $tmp{'amber'} ? $tmp{'amber'} : '',
	'TAA' => defined $tmp{'ochre'} ? $tmp{'ochre'} : '',
	'TGA' => defined $tmp{'opal'} ? $tmp{'opal'} : '',
);


# Meat of the tool
use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;
my $seqio_object = $bio->getSeqIO($options->{file});
my @fixed_seqs;
my @mod_seqs;
while(my $seq_object = $seqio_object->next_seq){
	foreach my $feat($seq_object->get_SeqFeatures){
		if($feat->primary_tag eq 'CDS'){
			# If we should ignore the stop that's currently in use
			my $stop = substr($bio->get_seq_from_feature($feat), -3);
			if($stop_ignore_map{$stop}){
				printf "Feature %s identified as a target\n", $bio->_getIdentifier($feat);
				# We want to get the next N bases
				my $extend_dist = 100;
				# Var to store if we have/not found a new stop
				my $new_end = -1;
				while($new_end == -1){
					# Get downstream
					my $downstream = $bio->get_seq_from_feature($feat,
						downstream => $extend_dist,
						parent => $seq_object,
					);
					# Remove original feature sequence
					$downstream = substr($downstream, $feat->end - $feat->start + 1);
					# Loop across stops
					for(my $i = 0; $i < length($downstream); $i += 3){
						my $codon = substr($downstream, $i, 3);
						# Check for valid stop
						if(defined $stop_ignore_map{$codon} && !$stop_ignore_map{$codon}){
							$new_end = $i;
							print "$codon ($new_end) ";
							$i = length($downstream);
						}
					}
					# If we haven't found an end, double the extend distance and check that much more sequence
					if($new_end == -1){
						$extend_dist *= 2;
					}
					# To prevent craziness, don't go overboard with how much downstream we'll check.
					if($extend_dist > 100_000){
						print STDERR $bio->_getIdentifier($feat) . " went out of range!!\n";
						last;
					}
				}
				# If we've found a new end, extend the feature (and include the stop codon)
				if($new_end > -1){
					my $location_type = ref $feat->location;


					if ( $location_type eq 'Bio::Location::Simple' ) {
						if($feat->strand() > 0){
							$feat->end($feat->end() + $new_end + 3);
						}else{
							$feat->start($feat->start() - $new_end - 3);
						}
					}
					elsif ( $location_type eq 'Bio::Location::Split' ) {
						my @sublocs = $feat->location->each_Location();
						foreach my $loc(@sublocs){
							if($feat->strand() > 0){
								$loc->end($loc->end() + $new_end + 3);
							}else{
								$loc->start($loc->start() - $new_end - 3);
							}
						}

					}else{
						warn "I don't know how to handle $location_type types of locations. Please submit a bug with Eric Rasche";
					}

					print "Extened " . $bio->_getIdentifier($feat) . " by  " . $new_end. "\n";
					my $custom_seq = $bio->get_seq_from_feature($feat);
					
					if($options->{translate}){
						my @nts = unpack("(A3)*", $custom_seq);
						foreach(my $i = 0; $i < scalar(@nts); $i++){
							if(defined $stop_ignore_map{$nts[$i]} && $stop_ignore_map{$nts[$i]}){
								$nts[$i] = $stop_ignore_map{$nts[$i]};
							}else{
								$nts[$i] = $bio->translate($nts[$i]);
							}
						}
						$custom_seq = join('', @nts);
					}
					push(@mod_seqs,
						sprintf(">%s %s bases added\n%s\n", $bio->_getIdentifier($feat), ($new_end+3), $custom_seq)
					);
				}
			}
		}
	}
	push(@fixed_seqs, $seq_object);
}

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
	name => 'results',
	libCPT => $libCPT,
);
$output->CRR(data => \@fixed_seqs);

my $modseq = CPT::OutputFiles->new(
	name => 'modified_features',
	libCPT => $libCPT,
);
$modseq->CRR(data => join('',@mod_seqs));


=head1 NAME

Virtual Suppressor Strain

=head1 DESCRIPTION

Given (one or more, mutli-genome) genbank file(s), this tool will extend features past specific stop codon types, suppressing them. It will also produce a fasta DNA sequence for modified features, allowing further investigation of those features which changed.

=head1 Acknowledgements

This tool was developed with significant input from Dr. Jason Gill @ TAMU

=cut
