#!/usr/bin/perl
use strict;
use warnings;
use Bio::AlignIO;
use CPT::GalaxyGetOpt;
use Data::Dumper;

# Via https://metacpan.org/pod/Bio::AlignIO#Bio::AlignIO-new
my %formats = (
	"bl2seq"    => "Bl2seq Blast output",
	"clustalw"  => "clustalw (.aln) format",
	"emboss"    => "EMBOSS water and needle format",
	"fasta"     => "FASTA format",
	"maf"       => "Multiple Alignment Format",
	"mase"      => "mase (seaview) format",
	"mega"      => "MEGA format",
	"meme"      => "MEME format",
	"msf"       => "msf (GCG) format",
	"nexus"     => "Swofford et al NEXUS format",
	"pfam"      => "Pfam sequence alignment format",
	"phylip"    => "Felsenstein PHYLIP format",
	"prodom"    => "prodom (protein domain) format",
	"psi"       => "PSI-BLAST format",
	"selex"     => "selex (hmmer) format",
	"stockholm" => "stockholm format",
);

my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[ 'file', 'Input file',
			{
				validate => 'File/Input',
				required => 1,
			}
		],
		[ 'min_length', 'Number of AAs or Percentage (should end in %)',
			{
				validate => 'String',
				required => 1,
				default => '10%',
			}
		],
		['alignment_format', 'Format of the MSA', { validate => 'Option', default => 'clustalw', options=>\%formats, required => 1}],
	],
	'outputs' => [
		[
			'results',
			'Findings report',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'msa_insertion_results',
				data_format    => 'text/html',
				default_format => 'HTML',
			}
		],
		[
			'interesting_sequences',
			'Sequences',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'msa_insertion_report',
				data_format    => 'genomic/raw',
				default_format => 'Fasta',
			}
		],
	],
	'defaults' => [
		'appid'   => 'cpt.msa.analysis',
		'appname' => 'MSA Insertion Analysis',
		'appdesc' => 'locates inserted regions in MSAs where a small number of sequences have large insertions',
		'appvers' => '1.94',
	],
	'tests' => [
		{
			test_name    => "Default",
			params => {
				'file' => 'test-data/inputs/pac-headful_1.msf',
				'alignment_format' => 'msf',
			},
			outputs => {
				'interesting_sequences' => ["msa_insertion_report.fa", 'test-data/outputs/msa_insertion_report.fa' ],
				'results' => ["msa_insertion_results.html", 'test-data/outputs/msa_insertion_results.html' ],
			},
		},
	],
);

my @files = glob("*.aln");
use CPT::Report::HTML;
my $report_html = CPT::Report::HTML->new();
my @report_fasta;

my $in  = Bio::AlignIO->new(-file   => $options->{file}, -format => $options->{alignment_format});
while ( my $aln = $in->next_aln() ) {
	if($options->{verbose}){
		printf STDERR "ALIGNMENT %s\n", $options->{file};
	}
	my $minimum_hit_length = calc_min($options->{min_length}, $aln->length());
	# Work backwards, by identifying regions that are commonly missed by sister sequences
	my %regions_of_interest;
	my @seqs;
	my $gap = $aln->gap_char();
	foreach my $seq ( $aln->each_seq() ){
		my @regions = identify_repeats($seq->seq(), $gap);
		foreach my $region(@regions){
			my ($start, $end) = @{$region};
			for(my $i = $start; $i <= $end; $i++){
				$regions_of_interest{$i}++;
			}
		}
		push(@seqs, $seq);
	}
	# Must be different from >50% of sisters, so we can zero out anything
	# less than 1/2 * num_seqs, and set everything else to 1. This should
	# help in identifying repeats.
	foreach(keys(%regions_of_interest)){
		if($regions_of_interest{$_} < int(scalar(@seqs)/2)){
			$regions_of_interest{$_} = 0;
		}else{
			$regions_of_interest{$_} = 1;
		}
	}
	# Create a region_str (111111000111101111100000)
	my $region_str = join('', map{ if(defined($regions_of_interest{$_})){$regions_of_interest{$_}}else{ 0 }} (0..$aln->length)); #keys(%regions_of_interest));
	# Find the repeated 1s, and filter (map) on minimum hit length
	my @region_locs = map { my ($s,$e) = @{$_}; if($e-$s > $minimum_hit_length){ $_ }else{ } } identify_repeats($region_str, '1+');
	# Loop across and identify which sequences are "SPECIAL".
	my $region_idx = 0;
	foreach(@region_locs){
		$region_idx++;
		my ($start, $end) = @{$_};
		#print "From $start to $end there is a large gap\n";
		# Looking for any sequences which has a positive # of residues in this location.
		my %count;
		foreach(@seqs){
			# Get the substr at this location
			my $substr = substr($_->seq, $start, $end-$start);
			# Remove all the gap characters, so we're left with coding chars
			$substr =~ s/$gap//g;
			# Count them
			$count{$_->display_id()}{num} = length($substr);
			$count{$_->display_id()}{seq} = substr($_->seq, $start, $end-$start);
			$count{$_->display_id()}{seq_strip} = $substr;
			$count{$_->display_id()}{region_len} = $end-$start;
		}
		# Print, sorted
		my $i = 0;
		my $html_alignment;

		foreach(sort{ $count{$b}{num} <=> $count{$a}{num} } keys(%count)){
			$html_alignment .= sprintf("%20s has %3d characters in this region [%s]\n", $_, $count{$_}{num}, $count{$_}{seq});

			if($count{$_}{num} > .5 * $count{$_}{region_len}){
				push(@report_fasta,
					Bio::Seq->new(
						-display_id => sprintf('%s_Region_%s_IDX_%s', $_, $region_idx, $i++),
						-seq => $count{$_}{seq_strip},
						-comment => sprintf('[%s AAs]', $count{$_}{num}),
					)
				);
			}
		}
		$report_html->a("<pre>" . $html_alignment . "</pre>");
	}
}

use CPT::OutputFiles;
my $out_html  = CPT::OutputFiles->new( GGO => $ggo, name => 'results' );
my $out_fasta = CPT::OutputFiles->new( GGO => $ggo, name => 'interesting_sequences' );

$out_html->CRR(data => $report_html->get_content());
$out_fasta->CRR(data => \@report_fasta);

sub identify_repeats {
	my ($seq, $gap) = @_;
	my @ret;
	while($seq =~ /$gap/g){
		push(@ret, [ $-[0], $+[0] ]);
	}
	return @ret;
}

sub calc_min {
	my ($string, $length) = @_;
	if($string =~ /%$/){
		$string = substr($string,0,index($string, '%'));
		$string /= 100;
		return int($string * $length);
	}
	return $string;
}

=head1 DESCRIPTION

This tool tries to find inserted sequences in an MSA. Insertions are limited to stretches longer that C<min_length> amino acids (or %), that does not align with the majority (>50%) of the others in the MSA.

=cut
