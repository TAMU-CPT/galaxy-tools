#!/usr/bin/perl

use strict;
use warnings;

use CPT::GalaxyGetOpt;
use Data::Dumper;

my $ggo  = CPT::GalaxyGetOpt->new();
my $options = $ggo->getOptions(
	'options' => [
		[ 'file|f', 'Input file',
			{
				validate => 'File/Input',
				#file_format => ['Genbank'],
				required => 1,
				file_format => ['genbank', 'embl', 'txt'],
			}
		],
	],
	'outputs' => [
		[
			'results',
			'Genbank Info and Statistics',
			{
				validate       => 'File/Output',
				required       => 1,
				default        => 'info',
				data_format    => 'text/html',
				default_format => 'HTML',
			}
		],
	],
	'defaults' => [
		'appid'   => 'PhageInfo',
		'appname' => 'Genbank Genome Information',
		'appdesc' => 'calculates various statistics and information about a genbank file',
		'appvers' => '1.94',
	],
	'tests' => [
	],
);

use CPT::Report::HTML;
my $report = CPT::Report::HTML->new();

use CPT::Bio;
my $bio = CPT::Bio->new();
use Bio::SeqIO;




my %wanted_tags = (CDS => 1);
my %pigs_data = {};


my $seqio_object2 = $bio->getSeqIO($options->{file});
while(my $seq_object2 = $seqio_object2->next_seq){
	# Coverage history, we'll count 0s vs 1s after
	my @coverage;

	# For all of our features, if we want them, bump coverage #s
	for my $feat_object ($seq_object2->get_SeqFeatures) {
		if($wanted_tags{ $feat_object->primary_tag }){
			my $loc = $feat_object->location;
			if(ref($loc) eq 'Bio::Location::Simple'){
				for(my $i=$feat_object->start; $i < $feat_object->end; $i++){
					$coverage[$i]++;
				}
			}elsif(ref($loc) ne 'Bio::Location::Fuzzy'){
				for my $location ( $loc->sub_Location ) {
					for(my $i=$location->start; $i < $location->end; $i++){
						$coverage[$i]++;
					}
				}
			}
		}
	}

	# Given that data, calculate coverages.
	my ($start, $end) = (1, $seq_object2->length());
	my ($pigs_covered, $art_covered) = (0,0);
	for(my $i=$start;$i<$end;$i++){
		if(defined $coverage[$i] && $coverage[$i] > 0){
			$pigs_covered++;
			$art_covered += $coverage[$i];
		}
	}
    $pigs_data{$seq_object2->display_id()} = [100*$art_covered/$end, 100*$pigs_covered/$end];
}

my $seqio_object = $bio->getSeqIO($options->{file});
# For all genomes in the GBK file
while(my $seq_object = $seqio_object->next_seq){
	my @features = $seq_object->get_SeqFeatures();
	my @genes;
	my $bases = 0;
	my @density_map;
	use Statistics::Descriptive;
	my $stat = Statistics::Descriptive::Full->new();
	my %g;
    my @coverage_results;
    my @coverage;
	foreach(@features){
		if($_->primary_tag eq 'CDS'){
			my %tags = map {$_ => 1} $_->get_all_tags();
			if(!defined($tags{'pseudo'}) && !defined($tags{'pseudogene'})){
				push(@genes, $_);
				my $len = $_->end() - $_->start();
				$bases += $len;
				$stat->add_data($len);
				for(my $i = int($_->start/1000); $i <= int($_->end/1000); $i++){
					$density_map[$i] += 1;
				}
				foreach my $b(split //,$_->seq->seq){
					$g{$b}++;
				}
			}
		}


		if($wanted_tags{ $_->primary_tag }){
			my $loc = $_->location;
			if(ref($loc) eq 'Bio::Location::Simple'){
				for(my $i=$_->start; $i < $_->end; $i++){
					$coverage[$i]++;
				}
			}elsif(ref($loc) ne 'Bio::Location::Fuzzy'){
				for my $location ( $loc->sub_Location ) {
					for(my $i=$location->start; $i < $location->end; $i++){
						$coverage[$i]++;
					}
				}
			}
		}



	}



    my ($start, $end) = (1, $seq_object->length());
    my ($pigs_covered, $art_covered) = (0,0);
    for(my $i=$start;$i<$end;$i++){
        if(defined $coverage[$i] && $coverage[$i] > 0){
            $pigs_covered++;
            $art_covered += $coverage[$i];
        }
    }
    my $cov_a = 100*$art_covered/$end;
    my $cov_b = 100*$pigs_covered/$end;



	my $density = 0;
	foreach(my $i=0;$i<scalar@density_map;$i++){
		if(defined $density_map[$i]){
			$density += $density_map[$i];
		}
	}
	$density /= scalar @density_map;
	my %go;
	foreach my $b(split //, $seq_object->seq){
		$go{$b}++;
	}

	my $g_sum = 0;
	foreach(keys(%g)){
		$g_sum+=$g{$_};
	}
	my $go_sum = 0;
	foreach(keys(%go)){
		$go_sum+=$go{$_};
	}








	$report->h1("Overview of " . $seq_object->display_id());
	$report->list_start('bullet');
	$report->list_element('Number of bases: ' . $seq_object->length());
	$report->list_element('Number of features: ' . scalar(@features));
	$report->list_end();

	$report->h3("Genes (CDS features withouth a /pseudo or /pseudogene qualifier)");
	$report->list_start('bullet');
	$report->list_element('count: ' . scalar(@genes));
	$report->list_element('bases: ' . $bases);
	$report->list_element(sprintf('density %0.3f genes per kb', $density));
	$report->list_element('average length: ' . int($stat->mean));
	$report->list_element('gene sequence composition:');
		$report->list_start('bullet');
		foreach(sort(keys(%g))){
			$report->list_element(sprintf("%s content: %s (%0.3f%%)", $_ , $g{$_} , 100*$g{$_}/$g_sum ));
		}
		$report->list_end();
	$report->list_element(sprintf('GC Percentage: %0.3f', 100*($g{'G'} + $g{'C'})/$g_sum ));
	$report->list_end();

	$report->h3('Overall sequence composition');
	$report->list_start('bullet');
	foreach(sort(keys(%g))){
		$report->list_element(sprintf("%s content: %s (%0.3f%%)", $_ , $go{$_} , 100*$go{$_}/$go_sum ));
	}
	$report->list_element(sprintf('GC Percentage: %0.3f', 100*($go{'G'} + $go{'C'})/$go_sum ));
	$report->list_end();

	$report->h3("Coverage");
	$report->list_start('bullet');
	$report->list_element('Genomic Coverage (CPT, not Artemis) ' . $cov_b);
	$report->list_end();

}

use CPT::OutputFiles;
my $output = CPT::OutputFiles->new(
        name => 'results',
        GGO => $ggo,
);
$output->CRR(data => $report->get_content());
