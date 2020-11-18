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
use Bio::Tools::SeqStats;
use CPT::BioData;
use CPT::Bio;

my $libCPT = CPT->new();

my $options = $libCPT->getOptions(
	'options' => [
		[
			'file' => 'Input Genbank file',
			{
				required => 1,
				validate => 'File/Input',
				multiple => 1,
			}
		],
	],
	'outputs' => [
		[
			'results',
			'Analysis results',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'trna',
				data_format    => 'text/tabular',
				default_format => 'CSV'
			}
		],
	],
	'defaults' => [
		'appid'   => 'tRNA',
		'appname' => 'tRNA Usage',
		'appdesc' => 'attempts to determine how many tRNAs are in the genome for a given amino acid',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);


# Preparing table header
my @header = ('Sequence');

# Used in table header
my $bd    = CPT::BioData->new();
my $bio   = CPT::Bio->new();
my %table = %{$bd->getTranslationTable()};

my %short_long_lookup = (
	'G',  => 'Gly',
	'P', => 'Pro',
	'A', => 'Ala',
	'V', => 'Val',
	'L', => 'Leu',
	'I', => 'Ile',
	'M', => 'Met',
	'C', => 'Cys',
	'F', => 'Phe',
	'Y', => 'Tyr',
	'W', => 'Trp',
	'H', => 'His',
	'K', => 'Lys',
	'R', => 'Arg',
	'Q', => 'Gln',
	'N', => 'Asn',
	'E', => 'Glu',
	'D', => 'Asp',
	'S', => 'Ser',
	'T', => 'Thr',
	'X', => 'XXX',
	'*' => 'Stop',
);


# Keys for header and lookup
my %aas;
foreach ( sort(keys(%table)) ) {
	$aas{$table{$_}} = 1;
} # gets XXX => X

my @gk = sort(keys(%aas));



foreach(sort(keys(%aas))){
	if($short_long_lookup{$_}){
		push(@header, sprintf('%s (%s)', $short_long_lookup{$_}, $_));
	}else{
		push(@header, $_);
	}
}

# https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/api/asn2gnb4.c
# Diff of that with my list of codons produced the following "extras"
foreach my $i('Asx','Gap','Glx','OTHER','Pro','Pyl','Sec','TERM','Xle'){
	$table{$i} = $i;
	push(@gk, $i);
	push(@header, $i);
}

my %global_results;
use Bio::SeqIO;
my @codon_usage_data;
foreach my $file(@{$options->{file}}){
	# For all genomes in the GBK file
	my $seqio_object = Bio::SeqIO->new(-file => $file,-format=>'genbank');

	while(my $seqobj = $seqio_object->next_seq()){
		$global_results{$seqobj->display_id()} = {};
		foreach my $feat($seqobj->get_SeqFeatures()){
			if($feat->primary_tag() eq 'tRNA'){
				my $codon_recognised;
				if($feat->has_tag('note')){
					my @vals = $feat->get_tag_values('note');
					foreach(@vals){
						if($_ =~ /codon_recognized=([ACTG]{3})/){
							$codon_recognised = $1;
						}
					}
				}
				#tRNA-Cys
				if($feat->has_tag('product')){
					my @vals = $feat->get_tag_values('product');
					foreach(@vals){
						if($_ =~ /tRNA-([A-Za-z]*)/){
							$codon_recognised = $bd->decode321($1);
							#Check with aa table. Maybe single letter code (tRNA-G) or one of the special ones (tRNA-Asx)
							if(!defined($codon_recognised)){
								$codon_recognised = $1;
							}
						}
					}
				}
				if(!defined($codon_recognised)){
					my $tag_list;
					foreach my $t($feat->get_all_tags()){
						foreach my $v($feat->get_tag_values($t)){
							$tag_list .= sprintf("\t/%s = '%s'\n", $t, $v);
						}
					}
					printf STDERR "Could not recognise a corresponding codon for tRNA [%s, %s, %s] in %s\n%s\n\n", $feat->start, $feat->end, $feat->strand, $seqobj->display_id(), $tag_list;
				}else{
					$global_results{$seqobj->display_id()}{$codon_recognised}++;
				}
				#print $feat,"\t$codon_recognised\n";
			}
		}
	}
}

my @data;

foreach my $seqid(sort(keys(%global_results))){
	my @row = ($seqid);
	my %local_trnas = %{$global_results{$seqid}};
	foreach my $trna(@gk){
		# If that start codon is used in that genome
		if(defined $local_trnas{$trna}){
			push(@row, $local_trnas{$trna});
		}else{
			push(@row, 0);
		}
	}
	push(@data, \@row);
}

my %results = (
	'Sheet1' => {
		header => \@header,
		data => \@data,
	}
);

use CPT::OutputFiles;
my $csv_output = CPT::OutputFiles->new(
       name => 'results',
       libCPT => $libCPT,
);
$csv_output->CRR(data => \%results);


=head1 NAME

tRNA usage

=head1 DESCRIPTION

Attempts to parse tRNAs listed in the genome and obtain number of codons for each amino acid.

=cut
